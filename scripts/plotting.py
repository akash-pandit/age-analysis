#!/usr/bin/env python

from scanpy.plotting import embedding
from matplotlib import pyplot as plt

from anndata import AnnData
from pathlib import Path


def plot_dual_celltypes(
        figtitle: str, 
        adata1: AnnData, 
        adata2: AnnData, 
        color_col1: str, 
        color_col2: str, 
        title1: str | None = None,
        title2: str | None = None,
        basis1: str = 'umap',
        basis2: str = 'umap',
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
        adata=adata1, 
        basis=basis1, 
        color=color_col1, 
        ax=axs[0], 
        show=False
    )
    axs[1] = embedding(
        adata=adata2, 
        basis=basis2,
        color=color_col2, 
        ax=axs[1], 
        show=False, 
        legend_loc=legend_loc
    )

    title1 = title1 or color_col1  # defaults to color_col1 if title1 is None
    title2 = title2 or color_col2

    axs[0].set_title(title1)
    axs[1].set_title(title2)

    if savepath:
        fig.savefig(Path(savepath))
