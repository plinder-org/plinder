# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import logging
from argparse import ArgumentParser
from os import environ
from pathlib import Path
from shlex import split
from shutil import which
from subprocess import PIPE, check_output
from sys import argv as sys_argv
from typing import List, Optional

from proc import Proc

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger("doc")
DEPENDENCY_BLOCKS = "test,pipeline,plots"


def get_dev_tag() -> str:
    """
    Get the current git describe which includes base image tag
    and whether or not the repository is in a dirty (uncommitted)
    state.

    Returns
    -------
    tag : str
        current git describe
    """
    return check_output(
        [
            "git",
            "describe",
            "--tags",
            "--always",
            "--dirty",
            "--abbrev=8",
        ],
        text=True,
    ).strip()


def get_version_bump(base_tag: str | None = None) -> str:
    """
    Inspect the git history for "bumpversion {phrase}"
    and return the new version.
    """
    if base_tag is None:
        current = get_dev_tag()
        base_tag = current.split("-")[0]
    bump = "patch"
    tokens = [
        "bumpversion major",
        "bumpversion minor",
        "bumpversion skip",
    ]
    for token in tokens:
        log = check_output(
            split(
                " ".join(
                    [
                        "git",
                        "log",
                        f"{base_tag}..HEAD",
                        "--oneline",
                        "--grep",
                        f"'{token}'",
                        "--format=%s",
                    ]
                )
            ),
            text=True,
        ).strip()
        # origin/main hack was to get around docs branch but isn'try:
        # present before pushing the tag back to main so just search for PR number
        if log and "(#" in log:
            bump = token.split()[1]
            break
    if bump == "skip":
        return ""
    try:
        import semver

        new_version = getattr(semver, f"bump_{bump}")(base_tag.lstrip("v"))
        new_version = f"v{new_version}"
        print(new_version)
        return new_version
    except ImportError:
        LOG.error("could not import semver")
        return ""


def get_image() -> str:
    """
    Minimize the diff if any of this changes
    """
    url = "ghcr.io"
    organization = "plinder-org"
    image_name = "plinder"
    return f"{url}/{organization}/{image_name}"


def get_docker() -> str:
    """
    Sanity check for docker in environment
    """
    docker = which("docker")
    if docker is None:
        raise ValueError("could not find docker executable!")
    if not (Path.cwd() / "docker-compose.yml").is_file():
        raise ValueError("could not find docker-compose.yml, are you in repo root?")
    return docker


def get_env(tag: str | None = None) -> dict[str, str]:
    """
    Get the build env that simulates docker compose

    Parameters
    ----------
    tag : str, default=None
        the tag to use instead of the current dev tag

    Returns
    -------
    env : dict[str, str]
        the build environment
    """
    image_repo = "/".join(get_image().split("/")[:-1])
    dev_tag = get_dev_tag()
    base_tag = dev_tag.split("-")[0]
    build_tag = tag or dev_tag
    LOG.info(f"get_env: base_tag={base_tag} build_tag={build_tag}")
    return {
        **environ,
        **{
            "BASE_TAG": base_tag,
            "IMAGE_REPO": image_repo,
            "BUILD_TAG": build_tag,
            "INDEX_URL": "https://pypi.org/simple",
        },
    }


def pull_base_image(build: bool = False, promote: bool = False) -> None:
    """
    Make sure the base image is up to date. Optionally
    promote the image tag using semver by inspecting the
    git history for "bumpverson {phrase}".

    Parameters
    ----------
    promote : bool, default=False
        if True, promote the image tag using semver
    """
    docker = get_docker()
    env = get_env()
    cmd = [docker, "compose", "pull", "base", "--quiet"]
    Proc(cmd, env=env).execute()
    if build:
        cmd = [
            docker,
            "build",
            "-f",
            "dockerfiles/base/Dockerfile",
            "-t",
            f"{env['IMAGE_REPO']}/plinder-base:{env['BASE_TAG']}",
            "--secret",
            "id=INDEX_URL",
            "--build-arg",
            f"BASE_TAG={env['BASE_TAG']}",
            "--build-arg",
            f"DEPENDENCY_BLOCKS={DEPENDENCY_BLOCKS}",
            (Path.cwd() / "dockerfiles/base").as_posix(),
        ]
        Proc(cmd, env=env).execute()
    if promote:
        promote_base_image(env)


def promote_base_image(env: dict[str, str]) -> None:
    """
    Promote the base image
    """
    base_tag = env["BASE_TAG"]
    new_tag = get_version_bump(base_tag)
    if not new_tag:
        LOG.info("skipping base image promotion")
        return
    image = f"{env['IMAGE_REPO']}/plinder-base"
    orig_image = f"{image}:{base_tag}"
    new_image = f"{image}:{new_tag}"
    LOG.info(f"promoting base image to {new_tag}")
    docker = get_docker()
    Proc([docker, "tag", orig_image, new_image], env=env).execute()
    check_output(["git", "tag", new_tag])


def build_image(tag: str | None = None, push: bool = False) -> str:
    """
    Build the image

    Parameters
    ----------
    tag : str, default=None
        the tag to use instead of the current dev tag
    push : bool, default=False
        push the image to the artifact registry

    Returns
    -------
    image : str
        the full image including tag
    """
    check_output(["git", "fetch", "--tags"])
    docker = get_docker()
    env = get_env(tag)
    image = f"{env['IMAGE_REPO']}/plinder"
    build_tag = env["BUILD_TAG"]
    registry = environ.get("PLINDER_REGISTRY", "/".join(get_image().split("/")[:-1]))
    cmd = [
        docker,
        "build",
        "-f",
        "dockerfiles/main/Dockerfile",
        "-t",
        f"{env['IMAGE_REPO']}/plinder:{build_tag}",
        "-t",
        f"{registry}/plinder:{build_tag}",
        "--secret",
        "id=INDEX_URL",
        "--build-arg",
        f"BASE_IMAGE={env['IMAGE_REPO']}/plinder-base",
        "--build-arg",
        f"BASE_TAG={env['BASE_TAG']}",
        "--build-arg",
        f"DEPENDENCY_BLOCKS={DEPENDENCY_BLOCKS}",
        Path.cwd().as_posix(),
    ]
    Proc(cmd, env=env).execute()
    if push:
        cmd = [docker, "push", f"{registry}/plinder:{build_tag}"]
        Proc(cmd, env=env).execute()
    return f"{image}:{build_tag}"


def test_image(
    tag: str,
    push: bool = False,
    args: Optional[List[str]] = None,
    dirty: bool = False,
) -> None:
    """
    Run the test service from docker compose. Optionally
    push both the base and plinder images to the artifact
    registry.

    Parameters
    ----------
    tag : str
        the tag to use instead of the current dev tag
    push : bool, default=False
        if True, push images to the artifact registry
    args : List[str], default=None
        the arguments to pass to the image
    dirty : bool, default=False
        if True, mount the current working tree
    """
    env = get_env(tag)
    docker = get_docker()
    cmd = [docker, "compose", "run", "-e", "PLINDER_OFFLINE=true"]
    if dirty:
        cmd.extend(["-v", f"{Path.cwd()}/src/plinder:/opt/conda/lib/python3.9/site-packages/plinder"])
    cmd.append("test")
    if args is not None and len(args):
        cmd.extend(
            split(f'''/bin/bash -c "python -m pytest -n auto -v {' '.join(args)} && cp .coverage reports/.coverage"''')
        )
    Proc(cmd, env=env).execute()
    if push:
        cmd = [docker, "compose", "push", "base", "plinder", "--quiet"]
        Proc(cmd, env=env).execute()


def run_image(
    args: Optional[List[str]] = None,
    build: bool = False,
    tag: str | None = None,
    dirty: bool = False,
    it: bool = False,
    script: str | None = None,
) -> None:
    """
    Run the image mounting the current working tree. This can be
    convenient for debugging during development.

    Parameters
    ----------
    image : str
        the image to run
    args : List[str], default=None
        the arguments to pass to the image
    """
    tag = tag or "latest"
    image = build_image(tag=tag) if build else get_image() + ":" + tag
    docker = which("docker")
    if docker is None:
        raise ValueError("could not find docker executable!")
    home = Path.home()
    host = Path.cwd() / "src"
    guest = "/opt/conda/lib/python3.10/site-packages/plinder"
    app = "/app/src"
    if not host.is_dir():
        raise RuntimeError(f"could not find {host}, please run from repo root")
    cmd = [
        docker,
        "run",
        "--platform",
        "linux/amd64",
        "-e",
        "PLINDER_TEST_BASE_DIR=/app",
        "-v",
        f"{host.parent}/tests:/app/tests/",
        "-v",
        f"{host.parent}/column_descriptions:/app/column_descriptions/",
        "-e",
        "PLINDER_MOUNT=/plinder",
        "-e",
        "PLINDER_LOG_LEVEL=10",
        "-v",
        f"{home}/.local/share/plinder:/plinder",
        "-v",
        f"{host}/plinder:{guest}",
        "-v",
        f"{host}:{app}",
    ]
    if script:
        cmd.extend(["-v", f"{host.parent}/{script}:{Path(app).parent}/{script}"])
    if it:
        import pty

        pty.spawn(cmd + ["-it", image, "/bin/bash"])
        return
    cmd.append(image)
    if args is not None:
        cmd.extend(args)
    Proc(cmd).execute()


def main(argv: Optional[List[str]] = None):
    """
    The main entrypoint for the docker helper

    Parameters
    ----------
    argv : List[str], default=None
        the command line arguments
    """
    if argv is None:
        argv = sys_argv[1:]
    parser = ArgumentParser(prog="plinder docker helper")
    subs = parser.add_subparsers(help="subcommand", dest="command")
    subs.add_parser("bump", help="Get version bump")
    pull = subs.add_parser("pull", help="Pull the base image")
    pull.add_argument(
        "--promote",
        default=False,
        action="store_true",
        help="Promote the image tag using semver",
    )
    build = subs.add_parser("build", help="Build the app image")
    test = subs.add_parser("test", help="Test the app image")
    test.add_argument("--dirty", default=False, action="store_true", help="Mount current working tree")
    run = subs.add_parser("run", help="Run the app image")
    run.add_argument("--it", default=False, action="store_true", help="Run in interactive mode")
    run.add_argument("--script", default="", help="Script to run")
    run.add_argument("--dirty", default=False, action="store_true", help="Mount current working tree")
    for sub in [build, test, run]:
        sub.add_argument(
            "--tag", default=None, help="The image tag to pass to build_image",
        )
    for sub in [build, test]:
        sub.add_argument(
            "--push",
            default=False,
            action="store_true",
            help="Push the image to the artifact registry",
        )
    for sub in [pull, run]:
        sub.add_argument(
            "--build",
            default=False,
            action="store_true",
            help="Build the image before running",
        )
    ns, args = parser.parse_known_args()
    nsargs = vars(ns)
    if len(args):
        nsargs["args"] = args
    command = nsargs.pop("command")
    func = {
        "bump": get_version_bump,
        "pull": pull_base_image,
        "build": build_image,
        "test": test_image,
        "run": run_image,
    }
    kwargs = {} if command == "bump" else nsargs
    if command is None:
        parser.print_help()
        exit()
    func[command](**kwargs)


if __name__ == "__main__":
    main()
