# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
from __future__ import annotations

import inspect
import logging
import os

LOGGING_FORMAT: str = "%(asctime)s | %(name)s:%(lineno)d | %(levelname)s : %(message)s"
try:
    DEFAULT_LOGGING_LEVEL: int = int(os.getenv("PLINDER_LOG_LEVEL", "20"))
except ValueError:
    DEFAULT_LOGGING_LEVEL = logging.INFO


class PlinderLoggingError(Exception):
    pass


def setup_logger(
    logger_name: str | None = None,
    log_level: int = DEFAULT_LOGGING_LEVEL,
    log_file: str | None = None,
    propagate: bool = False,
) -> logging.Logger:
    """
    Setup logger for the module name as the logger name by default
    for easy tracing of what's happening in the code

    Parameters
    ----------
    logger_name : str
        Name of the logger
    log_level : int
        Log level
    log_file: str | None
        optional log file to write to
    propagate : bool
        propagate log events to parent loggers, default = False

    Returns
    -------
    logging.Logger:
        logger object

    Examples
    --------
    >>> logger = setup_logger("some_logger_name")
    >>> logger.name
    'some_logger_name'
    >>> logger.level
    20
    >>> logger = setup_logger(log_level=logging.DEBUG)
    >>> logger.name
    'log.py'
    >>> logger.level
    10
    """

    if logger_name is None:
        # Get module name as the logger name, this is copied from:
        # https://stackoverflow.com/questions/13699283/how-to-get-the-callers-filename-method-name-in-python
        frame = inspect.stack()[1]
        module = inspect.getmodule(frame[0])
        file_path = __file__ if module is None else module.__file__
        logger_name = os.path.basename(file_path) if file_path is not None else "log"

    # Set up logger with the given logger name
    logger = logging.getLogger(logger_name)
    # Check if logging level has been set externally otherwise first pass logger.level == 0 (NOTSET)
    set_level = not bool(logger.level)
    if set_level:
        logger.setLevel(log_level)
    handler = logging.StreamHandler()
    if set_level:
        handler.setLevel(log_level)
    formatter = logging.Formatter(LOGGING_FORMAT)
    handler.setFormatter(formatter)
    if not len(logger.handlers):
        logger.addHandler(handler)

    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        if set_level:
            file_handler.setLevel(log_level)
        if not [h for h in logger.handlers if h.__class__ == logging.FileHandler]:
            logger.addHandler(file_handler)
    logger.propagate = propagate

    return logger
