"""
Plotting utilities for wgd
Arthur Zwaenepoel - 2017
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd


def plot_selection(dists, output_file=None, alphas=None, ks_range=(0.1, 5), offset=5, title='Species genus', **kwargs):
    """
    Plot a panel from the `all.csv` output from the Ks analysis.

    :param dists: distribution(s), if multiple provide as list
    :param alphas: alpha values for the different distributions (will assign automatically if not provided)
    :param ks_range: Ks range to include for plotting
    :param offset: offset of axis
    :param title: panel title
    :param kwargs: keyword arguments for :py:func:`matplotlib.pyplot.hist`
    :return: :py:class:`matplotlib.pyplot.Figure` object
    """
    fig = plt.figure(figsize=(12, 12))

    if type(dists) != list:
        dists = [dists]
        alphas = [alphas]

    if not alphas or not alphas[0]:
        alphas = list(np.linspace(0.2, 1, len(dists)))

    for i in range(len(dists)):
        dists[i] = dists[i][dists[i]['Ks'] > ks_range[0]]
        dists[i] = dists[i][dists[i]['Ks'] < ks_range[1]]

    # ks
    ax = fig.add_subplot(221)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([dist['Ks'] for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(dists[i]['Ks'], bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$K_S$')

    # ka
    ax = fig.add_subplot(222)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([dist['Ka'] for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(dists[i]['Ka'], bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$K_A$')

    # log(ka )
    ax = fig.add_subplot(223)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log(dist['Ka']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log(dists[i]['Ka']), bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$ln(K_A)$')

    # log(w)
    ax = fig.add_subplot(224)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log(dist['Omega']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log(dists[i]['Omega']), bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$ln(\omega)$')
    fig.suptitle('${0}$ ${1}$'.format(title.split()[0], title.split()[1]))

    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches='tight')

    return fig