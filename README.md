# HSC/MKP Age Analysis

An analysis studying aged mouse hematopoietic stem cells and megakaryocyte progenitors in an attempt to identify reasons for age-induced non-canonical platelet generation.

- [HSC/MKP Age Analysis](#hscmkp-age-analysis)
    - [Background](#background)
    - [Objective](#objective)
  - [Prerequisites](#prerequisites)
  - [Setup](#setup)

### Background

Hematopoietic Stem Cells (HSCs) are the source of our blood cells in both humans and mice. With a mouse model, [the Forsberg lab at UC Santa Cruz has shown](https://doi.org/10.1016/j.cell.2024.04.018) that as a mouse grows old, their HSCs begin to differentiate into a different kind Megakaryocyte Progenitors (the cells that govern platelet production), denoted 'non-canonical MKPs' or 'ncMKPs', as opposed to normally-behaving or 'canonical' MKPs. Non-canonical create higher amounts of platelets, and those platelets are also more prone to clotting (described as 'enhanced thrombosis' in the literature).

![The graphical abstract for Poscablo 2024](/assets/poscablo24-graphical-abstract.jpg)

**The combination of the weakening that comes with age with more aggressive, more likely clotting paints a concerning picture for *any* organism that may experience this, and can lead to an increase in heart disease.**

### Objective

This analysis hypothesizes that there is a transcriptomic or proteomic cause for ncMKP production and aims to find candidates for that cause.

## Prerequisites

Basic BASH shell proficiency is assumed. Mac users, this is built-in to your terminal. Windows users, please 

1. **Terminal Access**: Mac users use "Terminal". Windows users, please [install WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) or another BASH-based terminal emulator of your choice (e.g. Git Bash). 

2. **Git**: used for version control, [install here](https://git-scm.com/install/).

3. **uv**: Our Python package manager. Ensures that code runs the same way, no matter the machine. [Install it here](https://docs.astral.sh/uv/getting-started/installation/).

## Setup

To clone (download) this repository's contents, navigate to your parent directory of choice and run:
```bash
git clone https://github.com/akash-pandit/age-analysis
cd age-analysis
```

To download all python dependencies and configure your environment, run:
```bash
uv sync  # yeah, its that easy.
```

**For Jupyter users:** Launch jupyter with `uv` to ensure it uses the correct environment and navigate to one of the given URLs:
```bash
uv run jupyter lab
```

<!-- **For VSCode Users**: VSCode sometimes struggles to recognize a UV environment. To force VSCode's internal jupyter server to recognize a uv python environment, run the following line:
```bash
uv run python -m ipykernel install --user --name hsc-mkp-age-analysis --display-name "Age Analysis"
``` -->
