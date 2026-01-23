# Project Template

A template for bioinformatics workflows for the Forsberg lab. Replace this with a short description of whatever analysis is going on here. Created by Akash Pandit.

- [Project Template](#project-template)
  - [Prerequisites](#prerequisites)
  - [Setup](#setup)

## Prerequisites

Basic BASH shell proficiency is assumed. Mac users, this is built-in to your terminal. Windows users, please 

1. **Terminal Access**: Mac users use "Terminal". Windows users, please [install WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) or another BASH-based terminal emulator of your choice (e.g. Git Bash). 

2. **Git**: used for version control, [install here](https://git-scm.com/install/).

3. **uv**: Our Python package manager. Ensures that code runs the same way, no matter the machine. [Install it here](https://docs.astral.sh/uv/getting-started/installation/).

## Setup

To clone (download) this repository's contents, navigate to your parent directory of choice and run:
```bash
git clone https://github.com/akash-pandit/bnfo-analysis-template
cd bnfo-analysis-template 
```

To download all python dependencies and configure your environment, run:
```bash
uv sync  # yeah, its that easy.
```

**For Jupyter users:** Launch jupyter with `uv` to ensure it uses the correct environment and navigate to one of the given URLs:
```bash
uv run jupyter lab
```

**For VSCode Users**: VSCode sometimes struggles to recognize a UV environment. To force VSCode's internal jupyter server to recognize a uv python environment, run the following line:
```bash
uv run python -m ipykernel install --user --name internal_env_name --display-name "Env Name"
```
