# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
#!/usr/bin/env python

import argparse
import logging
import multiprocessing
import os
import time
from itertools import repeat
from multiprocessing import cpu_count
from subprocess import check_output
from typing import List, Optional

logging.basicConfig(format="[MPQueue] %(levelname)s: %(message)s", level=logging.DEBUG)


class Task:
    """A task class"""

    def __init__(
        self, command: str, path: str = ".", timeout: Optional[int] = None
    ) -> None:
        self.command = command
        self.path = path
        if timeout is not None:
            timeout = int(timeout)
        self.timeout = timeout

    def run(self) -> None:
        """Runs a command in the given path"""
        try:
            check_output(self.command, shell=True, timeout=self.timeout, cwd=self.path)
        except Exception:
            pass


class MPQueue(object):
    """Distributes Tasks to Pool"""

    def __init__(
        self,
        tasks: List[Task],
        num_cpus: int = 0,
        maxtasksperchild: Optional[int] = None,
        chunksize: int = 1,
    ) -> None:
        try:
            self.num_processes = int(num_cpus)
            if self.num_processes < 1:
                raise ValueError()
        except (ValueError, TypeError):
            logging.warning(
                "Number of cores not specified or incorrect. Using all cores."
            )
            self.num_processes = cpu_count() - 1

        logging.info(f"MPQueue will use {self.num_processes} cores")

        self.tasks = tasks
        self.num_tasks = len(tasks)
        self.maxtasksperchild = maxtasksperchild
        if maxtasksperchild is not None:
            self.maxtasksperchild = int(maxtasksperchild)
        self.chunksize = chunksize

    def process(self) -> None:
        with multiprocessing.get_context("spawn").Pool(
            processes=self.num_processes,
            maxtasksperchild=self.maxtasksperchild,
        ) as pool:
            pool.map(
                MPQueue.run_task,
                zip(enumerate(self.tasks), repeat(self.num_tasks)),
                chunksize=self.chunksize,
            )

    @staticmethod
    def run_task(tup: tuple[tuple[int, Task], int]) -> None:
        (item, task), total = tup
        logging.info(f"running task {item}/{total}: {task.command}")
        t0 = time.time()
        task.run()
        t1 = time.time()
        logging.info(f"running task {item}/{total} took {(t1-t0):.2f}s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="mpqueue")
    parser.add_argument(
        "tasks_file_name",
        help="A file containing a task for each line",
        metavar="tasks_file_name",
    )
    parser.add_argument(
        "--cores",
        "-cores",
        "-c",
        help="CPU cores to use",
        dest="cores",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--maxtasksperchild",
        "-maxtask",
        "-t",
        help="Maximum tasks per child process",
        dest="maxtasksperchild",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--chunksize",
        "-chunk",
        "-s",
        help="Chunksize to pass to Pool.map",
        dest="chunksize",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--cwd",
        default=".",
    )
    parser.add_argument(
        "--timeout",
        help="Timeout in seconds (must be intable)",
        default=None,
    )
    args = parser.parse_args()
    with open(args.tasks_file_name) as handle:
        tasks = []
        for line in handle:
            if line and not line.startswith("#"):
                tasks.append(Task(line.rstrip(os.linesep), args.cwd, args.timeout))

        queue = MPQueue(tasks, args.cores, args.maxtasksperchild)
        queue.process()
