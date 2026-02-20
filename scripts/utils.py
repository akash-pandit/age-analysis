#!/usr/bin/env python

from scanpy.plotting import embedding
from matplotlib import pyplot as plt

from anndata import AnnData
from scanpy import pl
from pathlib import Path


def setup_validate_paths(required_paths: list[Path], output_dirs: list[Path]):
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


def plot_dual_celltypes(
    figtitle: str,
    adata1: AnnData,
    adata2: AnnData,
    color_col1: str,
    color_col2: str,
    title1: str | None = None,
    title2: str | None = None,
    basis1: str = "umap",
    basis2: str = "umap",
    savepath: str | Path | None = None,
):
    """
    Plot embeddings for two AnnData objects. Shared labels are condensed to one legend.

    :param figtitle: Overarching figure title
    :type figtitle: str
    :param adata1: Lefthand AnnData object to plot
    :type adata1: AnnData
    :param adata2: Righthand AnnData object to plot
    :type adata2: AnnData
    :param color_col1: Keys for `adata1` annotations of observations/cells or variables/genes, e.g., `'ann1'` or `['ann1', 'ann2']`.
    :type color_col1: str
    :param color_col2: `color_col1` for `adata2`
    :type color_col2: str
    :param title1: AnnData plot title, defaults to the annotation color
    :type title1: str (optional)
    :param title2: `title1` for `adata2`
    :type title2: str (optional)
    :param basis1: key in `adata.obsm` to use for plotting. Checks for `{basis}` and `X_{basis}`. Defaults to 'umap'.
    :type basis1: str (optional)
    :param basis2: `basis1` for `adata2`
    :type basis2: str (optional)
    """

    # if adata1 and adata2 share the same legend, only plot one of 'em
    legend_loc = (
        None
        if all(adata1.uns[color_col1 + "_colors"] == adata2.uns[color_col2 + "_colors"])
        else "right margin"
    )

    fig, axs = plt.subplots(nrows=1, ncols=2, constrained_layout=True, figsize=(14, 5))
    fig.suptitle(t=figtitle)

    axs[0] = embedding(
        adata=adata1, basis=basis1, color=color_col1, ax=axs[0], show=False
    )
    axs[1] = embedding(
        adata=adata2,
        basis=basis2,
        color=color_col2,
        ax=axs[1],
        show=False,
        legend_loc=legend_loc,
    )

    title1 = title1 or color_col1  # defaults to color_col1 if title1 is None
    title2 = title2 or color_col2

    axs[0].set_title(title1)
    axs[1].set_title(title2)

    if savepath:
        fig.savefig(Path(savepath))


def plot_batches(adata: AnnData, batch_col: str, suptitle: str = "Treatment") -> None:
    batches = adata.obs[batch_col].unique()
    fig, axs = plt.subplots(
        nrows=1,
        ncols=len(batches),
        figsize=(5 * len(batches), 5),
        constrained_layout=True,
    )
    fig.suptitle(suptitle, size="xx-large")
    for ax, batch in zip(axs, batches):
        pl.umap(
            adata=adata[adata.obs[batch_col] == batch],
            ax=ax,
            show=False,
            title=f"Batch: {batch}",
            legend_loc="best",
        )
