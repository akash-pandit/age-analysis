#!/usr/bin/env python

from pathlib import Path


def setup_validate_paths(
        required_paths: list[Path],
        output_dirs: list[Path]
    ):
    """
    Validates existance of paths in `required_paths` and creates directories in `output_dirs`

    :param required_paths: list of paths that are required for successful program execution; throws an error for any missing path
    :type required_paths: list[Path]
    :param output_dirs: list of directories to write output to; creates any not-already-created directory (and intermediates)  with default permissions
    :type output_dirs: list[Path]
    """
    
    for path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"Could not find required input: {str(path)}")

    for dir in output_dirs:
        if not dir.is_dir():
            dir.mkdir(parents=True, exist_ok=True)
            