import matplotlib as mpl
from matplotlib import pyplot as plt
#from cycler import cycler
import numpy as np

def get_chrom_grid(chrom_lens, parent_gridelement=None):
    """
    subplot_spec is the parent gird
    """
    if parent_gridelement is None:
        parent_gridelement = mpl.gridspec.GridSpec(1, 1)[0]
    gs = mpl.gridspec.GridSpecFromSubplotSpec(1, len(chrom_lens),
                                              width_ratios=chrom_lens.astype(float) / sum(chrom_lens),
                                              subplot_spec=parent_gridelement)
    return gs, parent_gridelement


def get_chrom_axes(grid=None, parent_gridelement=None, fig=None,
                   ylabel="", plot_xlabels=True, color=None, plotfun=None):
    """
    axes of the returned figure are the chromosomes
    plus one global axes across all chromosomes
    """
    if plotfun is None:
        plotfun = plt.plot
    if grid is None:
        grid, parent_gridelement = get_chrom_grid(parent_gridelement=parent_gridelement)
    if fig is None:
        fig = plt.figure(1, figsize=(20, 2))
    fig.subplots_adjust(wspace=0)
    # print color
    if color is None:
        print("changed to current color")
        ax = plt.gca()
        color = plt.rcParams["axes.prop_cycle"].by_key()["color"][0]
    color = converter.to_rgb(color)
    chroms = series.index.get_level_values(level=0).unique()
    for i, (chrom, gr) in enumerate(zip(chroms, grid)):
        if i == 0:
            ax = fig.add_subplot(gr)
        else:
            ax = fig.add_subplot(gr, sharey=ax)
            ax.set_yticks([])

        if i % 2 == 0:
            color1 = tuple(np.array(color) * 0.6)
            if plot_xlabels:
                ax.set_xlabel('Chr' + chrom.rsplit("_")[1])
        else:
            color1 = color
        # print color1
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_prop_cycler([color1])
    ylim = ax.get_ylim()
    # we could also get the parent element with gr.get_topmost_subplotspec()
    total_ax = fig.add_subplot(parent_gridelement)
    total_ax.set_frame_on(False)
    total_ax.spines['top'].set_visible(True)
    total_ax.set_ylim(ylim)
    total_ax.set_xticks([])
    if not ylabel and series.name is not None:
        ylabel = series.name
    total_ax.set_ylabel(ylabel, fontsize=14)
    if plot_xlabels:
        total_ax.set_xlabel("Chromosomes", labelpad=25, fontsize=14)
    return fig.get_axes()


converter = mpl.colors.ColorConverter()


def plot_chrom_series(series, chrom_len, plotfun=None, grid=None, parent_gridelement=None, fig=None,
                      ylabel="", plot_xlabels=True, color=None, title=None, rightlabel=None, **kwa):
    """

    chrom_len ... series or dict with chrom names as keys and chromosomes length as values
    axes of the returned figure are the chromosomes
    plus one global axes across all chromosomes


    """
    if plotfun is None:
        plotfun = plt.plot
    if grid is None:
        grid, parent_gridelement = get_chrom_grid(chrom_lens=chrom_len, parent_gridelement=parent_gridelement)
    if fig is None:
        fig = plt.figure(1, figsize=(20, 2))
    fig.subplots_adjust(wspace=0)
    # print color
    if color is None:
        ax = plt.gca()
        color = plt.rcParams["axes.prop_cycle"].by_key()["color"][0]
    color = converter.to_rgb(color)
    print(color)
    # chroms = series.index.get_level_values(level=0).unique()
    try:
        chroms = chrom_len.index
    except AttributeError:
        chroms = chrom_len.keys()

    axes = {}
    for i, (chrom, gr) in enumerate(zip(chroms, grid)):
        # print chrom
        if i == 0:
            ax = fig.add_subplot(gr)
        else:
            ax = fig.add_subplot(gr, sharey=ax)
            ax.set_yticks([])

        axes.update({chrom: ax})

        if i % 2 == 0:
            color1 = tuple(np.array(color) * 0.6)
            if plot_xlabels:
                ax.set_xlabel(chrom)
        else:
            color1 = color
        # print color1
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_xlim([0, chrom_len[chrom]])
        # if chrom == 'CAE3':
        # print series.ix[chrom].index.values[~np.isnan(series.ix[chrom].values)]
        try:
            chrom_series = series.loc[chrom]
        except KeyError:
            try:
                chrom_series = series.loc[str(chrom)]
            except KeyError:
                try:
                    chrom_series = series.loc[int(chrom)]
                except (KeyError, TypeError):
                    continue

        plotfun(chrom_series.index.values, chrom_series.values, '.', color=color1, rasterized=True, **kwa)
    ylim = ax.get_ylim()
    # we could also get the parent element with gr.get_topmost_subplotspec()
    # we don't share y here because otherwise the ticks are interdependent,
    # is there a way to turn tick on for only one of the shared axes?
    total_ax = fig.add_subplot(parent_gridelement)
    axes.update({"total_ax": total_ax})
    total_ax.set_frame_on(False)
    total_ax.spines['top'].set_visible(True)
    total_ax.set_ylim(ylim)
    total_ax.set_xticks([])
    if not ylabel and series.name is not None:
        ylabel = series.name
    total_ax.set_ylabel(ylabel, fontsize=14)
    total_ax.set_yticks

    if plot_xlabels:
        total_ax.set_xlabel("Chromosomes", labelpad=25, fontsize=14)

    if rightlabel is not None:
        ax2 = ax.twinx()
        axes.update({"right_ax": ax2})
        ax2.set_ylabel(rightlabel, color=color1, fontsize=16)
        ax2.set_yticks([])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)

    if title is not None:
        plt.title(title)  # , position=(0.5,1.02),va='bottom'

    return fig, axes, grid, parent_gridelement