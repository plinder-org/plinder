# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""Compute the next release tag from "bumpversion" markers in merged PR titles."""

from __future__ import annotations

from subprocess import check_output


def get_dev_tag() -> str:
    """Return the current ``git describe`` (latest tag, distance, dirty state)."""
    return check_output(
        ["git", "describe", "--tags", "--always", "--dirty", "--abbrev=8"],
        text=True,
    ).strip()


def get_version_bump(base_tag: str | None = None) -> str:
    """
    Inspect the git history since the last tag for "bumpversion {major,minor,skip}"
    in a merged PR title and return the new version, defaulting to a patch bump.
    An empty string means no release.
    """
    if base_tag is None:
        base_tag = get_dev_tag().split("-")[0]
    bump = "patch"
    for token in ["bumpversion major", "bumpversion minor", "bumpversion skip"]:
        log = check_output(
            [
                "git",
                "log",
                f"{base_tag}..HEAD",
                "--oneline",
                "--grep",
                token,
                "--format=%s",
            ],
            text=True,
        ).strip()
        # only count merged PRs, whose titles carry the PR number
        if log and "(#" in log:
            bump = token.split()[1]
            break
    if bump == "skip":
        return ""
    import semver

    return f"v{getattr(semver, f'bump_{bump}')(base_tag.lstrip('v'))}"


if __name__ == "__main__":
    print(get_version_bump())
