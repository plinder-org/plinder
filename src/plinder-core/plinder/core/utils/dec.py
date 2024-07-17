# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from functools import wraps
from time import time
from typing import Any, Callable, TypeVar

from plinder.core.utils.log import setup_logger

T = TypeVar("T")


def timeit(func: Callable[..., T]) -> Callable[..., T]:
    """
    Simple function timer decorator
    """

    @wraps(func)
    def wrapped(*args: Any, **kwargs: Any) -> T:
        log = setup_logger(".".join([func.__module__, func.__name__]))
        ts = time()
        result = None
        try:
            result = func(*args, **kwargs)
            log.info(f"runtime succeeded: {time() - ts:.2f}s")
        except Exception:
            log.error(f"runtime failed: {time() - ts:.2f}s")
            raise
        return result

    return wrapped
