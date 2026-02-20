#!/usr/bin/env python

from pathlib import Path
import pyrootutils

ROOT = pyrootutils.setup_root(Path.cwd(), indicator='.git')
ASSETS = ROOT / "assets"
DATA = ROOT / "data"
OUTPUT = ROOT / "output"
FIGURES = ROOT / "figures"
NOTEBOOKS = ROOT / "notebooks"
SCRIPTS = ROOT / "scripts"
