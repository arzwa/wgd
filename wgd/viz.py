# TODO: still needs to be updated
import plumbum as pb
import matplotlib

if not 'DISPLAY' in pb.local.env:
    matplotlib.use('Agg')  # use this backend when no X server
import matplotlib.pyplot as plt
import matplotlib
import logging
import numpy as np
import seaborn as sns
import pandas as pd


def node_averages(df):
    # note that this returns a df with fewer rows, i.e. one for every
    # node in the gene family trees.
    return df.groupby(["family", "node"])["dS"].mean()

def node_weights(df):
    # note that this returns a df with the same number of rows
    return 1 / df.groupby(["family", "node"])["dS"].transform('count')

def parse_filter(s):
    x = [x.strip() for x in s.split("<")]
    if len(x) != 3:
        raise(ValueError("invalid 'x < field < y' filter string"))
    return (x[1], float(x[0]), float(x[2]))

def parse_filters(filterstring):
    return [parse_filter(s) for s in filterstring.split(",")]

def apply_filters(df, filters):
    for key, lower, upper in filters:
        df = df[df[key] > lower]
        df = df[df[key] < upper]
    return df

_labels = {
        "dS" : "$K_\mathrm{S}$",
        "dN" : "$K_\mathrm{A}$",
        "dN/dS": "$\omega$"}

def default_plot(
        *args, 
        alphas=None,
        colors=None,
        weighted=True, 
        title="",
        ylabel="duplication events",
        **kwargs):
    """
    Make a figure of node-weighted histograms for multiple distributions and
    variables. Returns the figure object.
    
    !!! note: Assumes the data frames are filtered as desired. 
    """
    ndists = len(args)
    alphas = alphas or list(np.linspace(0.2, 1, ndists))
    colors = colors or ['black'] * ndists
    
    # assemble panels
    keys = ["dS", "dS", "dN", "dN/dS"]
    funs = [lambda x: x, np.log10, np.log10, np.log10]
    fig, axs = plt.subplots(2, 2)
    for (c, a, dist) in zip(colors, alphas, args):
        for ax, k, f in zip(axs.flatten(), keys, funs):
            w = node_weights(dist)
            x = f(dist[k])
            y = x[np.isfinite(x)]
            w = w[np.isfinite(x)]
            ax.hist(y, weights=w, color=c, alpha=a, **kwargs)
            xlabel = _labels[k]
            if f == np.log10:
                xlabel = "$\log_{10}" + xlabel[1:-1] + "$"
            ax.set_xlabel(xlabel)
    axs[0,0].set_ylabel(ylabel)
    axs[1,0].set_ylabel(ylabel)

    # finalize plot
    sns.despine(offset=1)
    fig.suptitle(title, x=0.1, y=0.99, ha="left", va="top")
    fig.tight_layout()
    plt.subplots_adjust(top=0.85)  # prevent suptitle from overlapping
    return fig


def syntenic_depth_plot(segprofile):
    cols = segprofile.columns
    n = len(cols)
    fig, axs = plt.subplots(1, int(n + n*(n-1)/2))
    if n == 1:
        axs = [axs]  # HACK
    k = 0
    for i in range(n):
        for j in range(i, n):
            pairs, counts = dupratios(segprofile[cols[i]], segprofile[cols[j]])
            ax = axs[k]
            ax.barh(np.arange(len(pairs)), counts, color="k", alpha=0.2)
            ax.set_yticks(np.arange(len(pairs)))
            ax.set_yticklabels(["{}:{}".format(int(x[0]), int(x[1])) for x in pairs])
            ax.set_title("{}:{}".format(cols[i], cols[j]), fontsize=9)
            k += 1
    for ax in axs:
        ymn, ymx = ax.get_ylim()
        ax.set_ylim(-0.5, ymx)
        ax.set_xlabel("# segments")
    axs[0].set_ylabel("A:B ratio")
    sns.despine(trim=False, offset=3)
    fig.tight_layout()
    return fig


def dupratios(col1, col2, by="first"):
    d = {}
    for pair in zip(col1,col2):
        if pair not in d:
            d[pair] = 0
        d[pair] += 1
    if by == "first":
        keyfun = lambda x: x
    elif by == "ratio":
        lambda x: x[0]/(1+x[1])
    elif by == "second":
        keyfun = lambda x: x[1]
    kys = sorted(d, key=keyfun)
    return kys, [d[k] for k in kys]
