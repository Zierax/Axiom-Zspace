"""
Logging configuration module for Axiom-ZSpace production pipeline.

This module provides centralized logging setup that reads configuration from
config/production.yaml and configures Python's logging module with appropriate
handlers for both console and file output.

Requirements: 6.4, 6.8
"""

import logging
import sys
import yaml
from pathlib import Path
from typing import Optional


class _UnicodeStreamHandler(logging.StreamHandler):
    """StreamHandler that gracefully handles Unicode on non-UTF-8 consoles.

    On Windows with legacy codepages (cp1256, cp1252, etc.) the default
    StreamHandler raises UnicodeEncodeError when log messages contain
    characters like U+2192 (arrow), U+0394 (delta), or U+00D7 (multiply).

    This handler overrides emit() to catch encoding errors and replace
    unencodable characters with '?' instead of crashing.
    """

    def emit(self, record):
        try:
            msg = self.format(record)
            stream = self.stream
            try:
                stream.write(msg + self.terminator)
            except UnicodeEncodeError:
                # Encode with 'replace' to handle non-encodable characters
                encoding = getattr(stream, 'encoding', 'utf-8') or 'utf-8'
                safe_msg = msg.encode(encoding, errors='replace').decode(encoding, errors='replace')
                stream.write(safe_msg + self.terminator)
            self.flush()
        except RecursionError:
            raise
        except Exception:
            self.handleError(record)


def suppress_astroquery_logger() -> None:
    """Silence astroquery logs WITHOUT creating its logger prematurely.

    Calling logging.getLogger('astroquery') before astroquery has initialized
    its own logger registers a plain logging.Logger under that name. When
    astroquery's _init_log() later runs, getLogger('astroquery') then returns
    that plain Logger instead of an AstropyLogger, and its ._set_defaults()
    call raises AttributeError. Only touch the logger if astroquery has
    already initialized it.
    """
    try:
        if 'astroquery' in logging.Logger.manager.loggerDict:
            logging.getLogger('astroquery').setLevel(logging.WARNING)
    except Exception:
        pass  # Silently ignore if astroquery logger setup fails


def setup_logging(config_path: str = "config/production.yaml") -> None:
    """
    Configure Python logging module based on production.yaml settings.
    
    Reads logging configuration from the specified YAML file and sets up
    logging handlers for console and/or file output. Supports configurable
    log levels (DEBUG, INFO, WARNING, ERROR).
    
    Args:
        config_path: Path to production configuration YAML file.
                    Defaults to "config/production.yaml".
    
    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is malformed
        ValueError: If log level is invalid
    
    Example:
        >>> from zspace_engine.logging_config import setup_logging
        >>> setup_logging()
        >>> import logging
        >>> logging.info("Pipeline started")
    """
    # Initialize astropy logging FIRST to ensure its custom logger class is set up
    # This prevents the "'Logger' object has no attribute '_set_defaults'" error
    try:
        from astropy import log as astropy_log
        # Just accessing it ensures astropy's logging is initialized
        _ = astropy_log.level
    except ImportError:
        pass  # astropy not installed, no problem
    
    # Load configuration
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
    except IOError as e:
        raise IOError(f"Failed to read configuration file {config_path}: {e}") from e
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Failed to parse YAML configuration file {config_path}: {e}") from e
    
    # Extract logging configuration
    log_config = config.get("logging", {})
    
    # Parse log level
    level_str = log_config.get("level", "INFO").upper()
    try:
        level = getattr(logging, level_str)
    except AttributeError:
        raise ValueError(
            f"Invalid log level: {level_str}. "
            f"Must be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL"
        )
    
    # Build handlers list
    handlers = []
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Console handler
    if log_config.get("console", True):
        console_handler = _UnicodeStreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        handlers.append(console_handler)
    
    # File handler
    log_file = log_config.get("file")
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Configure root logger WITHOUT removing existing handlers
    # This is critical to preserve astropy's custom logger class
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    
    # Add our handlers to root logger (don't remove existing ones)
    for h in handlers:
        # Check if similar handler already exists
        handler_exists = False
        for existing_h in root_logger.handlers:
            if type(existing_h) == type(h):
                if isinstance(h, logging.FileHandler):
                    if existing_h.baseFilename == h.baseFilename:
                        handler_exists = True
                        break
                elif isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler):
                    handler_exists = True
                    break
        
        if not handler_exists:
            root_logger.addHandler(h)
    
    # Suppress verbose astroquery logging (without creating its logger early)
    suppress_astroquery_logger()
    
    # Log successful configuration
    logger = logging.getLogger(__name__)
    logger.info(f"Logging configured: level={level_str}, console={log_config.get('console', True)}, file={log_file}")


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """
    Get a logger instance for a specific module.
    
    This is a convenience wrapper around logging.getLogger() that ensures
    setup_logging() has been called. If not, it will call it with defaults.
    
    Args:
        name: Logger name, typically __name__ of the calling module.
              If None, returns the root logger.
    
    Returns:
        Configured logger instance
    
    Example:
        >>> from zspace_engine.logging_config import get_logger
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing TIC 12345678")
    """
    # Check if logging has been configured
    root_logger = logging.getLogger()
    if not root_logger.handlers:
        # No handlers configured, set up with defaults
        try:
            setup_logging()
        except FileNotFoundError:
            # Fallback to basic config if production.yaml not found
            root = logging.getLogger()
            if not root.handlers:
                h = _UnicodeStreamHandler()
                h.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
                root.addHandler(h)
                root.setLevel(logging.INFO)
    
    # Disable verbose astroquery logging (without creating its logger early)
    suppress_astroquery_logger()
    
    return logging.getLogger(name)
