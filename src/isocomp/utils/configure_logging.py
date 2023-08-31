import logging
import logging.config
import os


def configure_logging(level=logging.INFO,
                      to_file: bool = False,
                      filename: str = None):
    """Configure logging to print to the console and/or a file.

    :param level: The logging level, defaults to logging.INFO
    :type level: int, optional
    :param to_file: Whether to log to a file, defaults to False
    :type to_file: bool, optional
    :param filename: The name of the file to log to, defaults to None
    :type filename: str, optional

    :raises ValueError: If the logging level is not valid
    :raises ValueError: If to_file is not a boolean
    :raises ValueError: If filename is not a string
    :raises ValueError: If filename is not provided when to_file is True
    :raises FileNotFoundError: If the directory for the log file does not exist

    Configure Logging Examples
    --------------------------
    >>> configure_logging()
    >>> configure_logging(level=logging.DEBUG)
    >>> configure_logging(level=logging.INFO)
    """
    if level not in [logging.DEBUG,
                     logging.INFO,
                     logging.WARNING,
                     logging.ERROR,
                     logging.CRITICAL]:
        raise ValueError("The logging level must be one of logging.DEBUG, "
                         "logging.INFO, logging.WARNING, logging.ERROR, "
                         "or logging.CRITICAL.")
    
    if not isinstance(to_file, bool):
        raise ValueError("to_file must be a boolean.")
    
    if to_file:
        if not filename:
            raise ValueError("filename must be provided when to_file is True.")
        if not isinstance(filename, str):
            raise ValueError("filename must be a string.")
        # the `or` operator returns the first value if it is truthy
        log_file_dir = os.path.dirname(filename) or '.'
        if not os.path.isdir(log_file_dir):
            raise FileNotFoundError("The directory for the log file does not "
                                    "exist.")

    handlers = {}

    if to_file and filename:
        handlers['file'] = {
            'class': 'logging.FileHandler',
            'filename': filename,
            'mode': 'a',
        }
    else:
        handlers['console'] = {
            'class': 'logging.StreamHandler',
        }

    LOGGING_CONFIG = {
        'version': 1,
        'disable_existing_loggers': False,
        'handlers': handlers,
        'root': {
            'handlers': list(handlers.keys()),
            'level': level,
        },
    }

    logging.config.dictConfig(LOGGING_CONFIG)
