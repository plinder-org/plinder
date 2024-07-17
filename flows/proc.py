# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
import io
import logging
from subprocess import PIPE, STDOUT, Popen
from typing import Any, Dict, List, Optional, Union

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger("proc")


class Proc(Popen):
    """
    An opinionated Popen with an execute method
    to mirror communicate that handles input as strings,
    supports early return, and logs as the process runs.
    Also retain stdout and stderr as a list of strings
    for downstream access.
    """

    def __init__(
        self,
        *args,
        stdin=PIPE,
        stdout=PIPE,
        stderr=STDOUT,
        text=True,
        **kwargs,
    ):
        super().__init__(
            *args,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            text=text,
            **kwargs,
        )
        self._stderr = stderr
        # squelch the static type checkers that
        # see these pipes as None since we're using
        # PIPEs for everything
        self.stdin: io.TextIOWrapper
        self.stdout: io.TextIOWrapper
        self.stderr: io.TextIOWrapper

    def execute(
        self,
        *,
        inputs: Optional[str] = None,
    ) -> None:
        """
        Like Popen.communicate but tail output

        Parameters
        ----------
        inputs : str
            a command to pass to the subprocess

        """
        if isinstance(self.args, list):
            LOG.info(" ".join(map(str, self.args)))
        else:
            LOG.info(self.args)

        streams = [("stdout", "stdout"), ("stderr", "stderr")]
        if self._stderr == STDOUT:
            streams = [("stdout", "stream")]
        with self:
            try:
                if inputs is not None:
                    self.stdin.write(inputs)
                    self.stdin.close()
                for attr, stream in streams:
                    for line in getattr(self, attr):
                        LOG.info(f"{stream}: {line.rstrip()}")
            except Exception:
                self.kill()
