#!/usr/bin/python3
"""
--------------------------------------------------------------------------------

Copyright (C) 2018 Arthur Zwaenepoel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact: arzwa@psb.vib-ugent.be

--------------------------------------------------------------------------------
"""
import numpy as np
import logging
from sklearn import mixture
import scipy.stats as ss
import plumbum as pb
import matplotlib

if not 'DISPLAY' in pb.local.env:
    matplotlib.use('Agg')  # use this backend when no X server
import matplotlib.pyplot as plt
import seaborn as sns


def filter_group_data(
        df, aln_id=0, aln_len=300, aln_cov=0, min_ks=0, max_ks=5,
        weights_outliers_included=False
):
    """
    Filter a Ks data frame for modling purposes.

    :param df: data frame
    :param aln_id: alignment identity threshold
    :param aln_len: alignment length threshold
    :param aln_cov: alignment coverage threshold
    :param min_ks: minimum Ks value
    :param max_ks: maximum Ks value
    :param weights_outliers_included: boolean, adapt weights when removing
        outliers
    :return: filtered data frame
    """
    # pre node-grouping filters, irrespective of outlier removal
    df = df.dropna()
    df = df[df["AlignmentCoverage"] >= aln_cov]
    df = df[df["AlignmentIdentity"] >= aln_id]
    df = df[df["AlignmentLength"] >= aln_len]

    # Ks range filters
    # if one filters before the node-weighting, the weights are adapted with
    # respect to the outlier filtering (as in Vanneste et al. 2013)
    if not weights_outliers_included:
        df = df[df["Ks"] > min_ks]
        df = df[df["Ks"] <= max_ks]

    # grouping
    df = df.groupby(['Family', 'Node']).mean()

    if weights_outliers_included:
        df = df[df["Ks"] > min_ks]
        df = df[df["Ks"] <= max_ks]

    return df


def get_array_for_mixture(df):
    """
    Get an array of log transformed Ks values.

    :param df: data frame
    :return: array
    """
    X = np.array(df["Ks"].dropna())
    X = X[X > 0]
    X = np.log(X).reshape(-1, 1)
    return X


def reflect(data):
    """
    Reflect an array around it's left boundary

    :param data: np.array
    :return: np.array
    """
    reflection = -1 * (data - np.min(data)) + np.min(data)
    return np.hstack([data, reflection])


def reflected_kde(df, min_ks, max_ks, bandwidth, bins, out_file):
    """
    Perform Kernel density estimation (KDE) with reflected data.
    The data frame is assumed to be grouped already by 'Node'.

    :param df: data frame
    :param min_ks: minimum Ks value (best is to use 0, for reflection purposes)
    :param max_ks: maximum Ks value
    :param bandwidth: bandwidth (None results in Scott's rule of thumb)
    :param bins: number of histogram bins
    :param out_file: output file
    :return: nada
    """
    ks = np.array(df['Ks'])
    ks_reflected = reflect(ks)
    fig, ax = plt.subplots(figsize=(9, 4))
    if bandwidth:
        ax = sns.distplot(
                ks_reflected, bins=bins * 2, ax=ax,
                hist_kws={"rwidth": 0.8, "color": "k", "alpha": 0.2},
                kde_kws={"bw": bandwidth},
                color="k"
        )
    else:
        ax = sns.distplot(
                ks_reflected, bins=bins * 2, ax=ax,
                hist_kws={"rwidth": 0.8, "color": "k", "alpha": 0.2},
                color="k"
        )
    ax.set_xlim(min_ks, max_ks)
    sns.despine(offset=5, trim=False)
    ax.set_ylabel("Density")
    ax.set_xlabel("$K_{\mathrm{S}}$")
    fig.savefig(out_file, bbox_inches='tight')


def fit_gmm(X, n1, n2, max_iter=100, n_init=1, **kwargs):
    """
    Compute Gaussian mixtures for different numbers of components

    :param X: data frame (log transformed Ks values)
    :param n1: minimum number of components
    :param n2: maximum number of components
    :param max_iter: maximum number of iterations
    :param n_init: number of k-means initializations
    :param kwargs: other keyword args for `GaussianMixture`
    :return: models, bic, aic, best model
    """
    # fit models with 1 to n components
    N = np.arange(n1, n2 + 1)
    models = [None for i in range(len(N))]
    for i in range(len(N)):
        logging.info("Fitting GMM with {} components".format(N[i]))
        models[i] = mixture.GaussianMixture(
                n_components=N[i], covariance_type='full', max_iter=max_iter,
                n_init=n_init, **kwargs
        ).fit(X)
        logging.info("Component mean, variance, weight: ")
        log_components(models[i])

    # compute the AIC and the BIC
    aic = [m.aic(X) for m in models]
    bic = [m.bic(X) for m in models]
    best = models[np.argmin(bic)]

    return models, bic, aic, best


def log_components(m):
    for j in range(len(m.means_)):
        logging.info(".. {0:.3f}, {1:.3f}, {2:.3f}".format(
                np.exp(m.means_[j][0]),
                m.covariances_[j][0][0], m.weights_[j])
        )


def fit_bgmm(X, n1, n2, gamma=1e-3, max_iter=100, n_init=1, **kwargs):
    """
    Compute Bayesian Gaussian mixture

    :param X: data frame (log transformed Ks values)
    :param n1: minimum number of components
    :param n2: maximum number of components
    :param gamma: inverse of regularization strength
    :param max_iter: maximum number of iterations
    :param n_init: number of k-means initializations
    :param kwargs: other keyword args for `GaussianMixture`
    :return: models
    """
    # fit models with 1 to n components
    N = np.arange(n1, n2 + 1)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        logging.info("Fitting BGMM with {} components".format(N[i]))
        models[i] = mixture.BayesianGaussianMixture(
                weight_concentration_prior=gamma, n_init=n_init,
                n_components=N[i], covariance_type='full', max_iter=max_iter,
                **kwargs
        ).fit(X)
        log_components(models[i])

    return models


def plot_mixture(model, data, ax, l=0, u=5, color='black', alpha=0.2,
                 log=False, bins=25, alpha_l1=1):
    """
    Plot a mixture model. Assumes a log-transformed model and data
    and will back-transform.

    Note from scipy docs:
    ---------------------
    A common parametrization for a lognormal random variable Y
    is in terms of the mean, mu, and standard deviation, sigma,
    of the unique normally distributed random variable X such
    that exp(X) = Y. This parametrization corresponds to setting
    s = sigma and scale = exp(mu).

    So since we fit a normal mixture on the logscale, we can get
    the lognormal by setting scale to np.exp(mean[k]) and s to
    np.sqrt(variance[k]) for component k.

    :param model: model object
    :param data: data array
    :param ax: figure ax
    :param l: lower Ks limit
    :param u: upper Ks limit
    :param color: color for histogram
    :param alpha: alpha value
    :param log: plot on log scale?
    :param bins: number of histogram bins
    :param alpha_l1: alpha value for mixture lines
    :return: ax
    """
    x = np.linspace(l, u, 1000).reshape((-1, 1))
    if not log:
        data = np.exp(data)
    data = data[data >= l]
    data = data[data <= u]
    bins = ax.hist(
            data, bins, rwidth=0.8, color=color, alpha=alpha, density=True)
    maxy = max(bins[0])
    means = model.means_
    varcs = model.covariances_
    weights = model.weights_
    mix = None
    first = True
    for k in range(len(means)):
        if not log:
            curve = ss.lognorm.pdf(
                    x, scale=np.exp(means[k]), s=np.sqrt(varcs[k])) * weights[k]
        else:
            curve = ss.norm.pdf(
                    x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k]
        ax.plot(x, curve, '--k', color='black', alpha=0.4)
        if first:
            mix = curve
            first = False
        else:
            mix += curve
    ax.plot(x, mix, color='black', alpha=alpha_l1)
    ax.set_xlim(l, u)
    if log:
        ax.set_xlabel("$\mathrm{log}(K_{\mathrm{S}})$")
    else:
        ax.set_xlabel("$K_{\mathrm{S}}$")
    ax.set_ylabel("Density")
    return ax


def inspect_aic(aic):
    """
    Evaluate the Akaike information criterion results for mixture models.

    :param aic: AIC values
    """
    im = np.argmin(aic)
    logging.info("")
    logging.info("AIC assessment:")
    logging.info("min(AIC) = {:.2f} for model {}".format(aic[im], im + 1))
    logging.info("Relative probabilities compared to model {}:".format(im + 1))
    logging.info("   /                          \\")
    logging.info("   |      (min(AIC) - AICi)/2 |")
    logging.info("   | p = e                    |")
    logging.info("   \                          /")
    for i, aic_i in enumerate(aic):
        p_i = np.exp((aic[im] - aic_i) / 2)
        logging.info(".. model{:>4}: p = {:.4f}".format(i + 1, p_i))
    logging.info("")


def inspect_bic(bic):
    """
    Evaluate the BIC values

    :param bic: BIC values
    """
    im = np.argmin(bic)
    l = [
        "0 to  2:   Very weak",
        "2 to  6:    Positive",
        "6 to 10:      Strong",
        "    >10: Very Strong"
    ]
    logging.info("")
    logging.info("Delta BIC assessment: ")
    logging.info("min(BIC) = {:.2f} for model {}".format(bic[im], im + 1))
    for i, bic_i in enumerate(bic):
        dbic = bic_i - bic[im]
        j = 0
        if dbic > 2: j = 1
        if dbic > 6: j = 2
        if dbic > 10: j = 3
        logging.info(".. model{:>4}: delta(BIC) = {:>8.2f} ({})".format(
                i + 1, dbic, l[j]))
    logging.info("")


def plot_probs(m, ax, l=0.0, u=5, ylab=True):
    """
    Plot posterior component probabilities

    :param m: model
    :param ax: figure ax
    :param l: lower Ks limit
    :param u: upper Ks limit
    :param ylab: plot ylabel?
    :return: ax
    """
    if l == 0:
        l = 0.005
    x = np.linspace(l, u, 1000).reshape(-1, 1)
    p = m.predict_proba(np.log(x))
    p = p.cumsum(1).T
    x = np.linspace(l, u, 1000)
    order = tuple(np.argsort(m.means_.reshape((1, -1)))[0])
    alphas = np.linspace(0.2, 1, p.shape[0])[order,]
    for i, array in enumerate(p):
        if i == 0:
            ax.fill_between(
                    x, 0, p[i], color='gray', alpha=alphas[i],
                    label=order[i] + 1)
        elif i == len(alphas) - 1:
            ax.fill_between(
                    x, p[i - 1], 1, color='gray', alpha=alphas[i],
                    label=order[i] + 1)
        else:
            ax.fill_between(
                    x, p[i - 1], p[i], color='gray', alpha=alphas[i],
                    label=order[i] + 1
            )
    ax.set_xlim(0, u)
    ax.set_ylim(0, 1)
    ax.set_xlabel('$K_{\mathrm{S}}$')
    if ylab:
        ax.set_ylabel('$P(class|K_{\mathrm{S}})$')
    ax.legend(frameon=True)
    return ax


def plot_aic_bic(aic, bic, min_n, max_n, out_file):
    """
    Plot AIC and BIC curves

    :param aic: aic values
    :param bic: bic values
    :param out_file: output file
    :return: nada
    """
    x_range = list(range(min_n, max_n + 1))
    fig, axes = plt.subplots(1, 2, figsize=(12, 3))
    axes[0].plot(np.arange(1, len(aic) + 1), aic, color='k', marker='o')
    axes[0].set_xticks(list(range(1, len(aic) + 1)))
    axes[0].set_xticklabels(x_range)
    axes[0].grid(ls=":")
    axes[0].set_ylabel("AIC")
    axes[0].set_xlabel("# components")
    axes[1].plot(np.arange(1, len(bic) + 1), bic, color='k', marker='o')
    axes[1].set_xticks(list(range(1, len(bic) + 1)))
    axes[1].set_xticklabels(x_range)
    axes[1].grid(ls=":")
    axes[1].set_ylabel("BIC")
    axes[1].set_xlabel("# components")
    fig.tight_layout()
    fig.savefig(out_file)


def plot_bars_weights(model, ax):
    """
    Plot Component weights

    :param model: model
    :param ax: figure ax
    :return: ax
    """
    ax.bar(np.arange(1, model.n_components + 1), model.weights_, color='k',
           alpha=0.2)
    for i, mn in enumerate(model.means_):
        ax.text(
                x=i + 1, y=model.weights_[i] + 0.1,
                horizontalalignment='center',
                s='$\hat{{\mu}} = {:.2f}$'.format(np.exp(mn)[0])
        )
    ax.set_ylabel("weight")
    ax.set_xlabel("component")
    ax.set_ylim(0, 1)
    return ax


def plot_all_models_gmm(models, data, l, u, bins, out_file):
    """
    Plot a bunch of GMMs.

    :param models: list of GMM model objects
    :param data: Ks array
    :param l: lower Ks limit
    :param u: upper Ks limit
    :param bins: number of histogram bins
    :param out_file: output file
    :return: nada
    """
    fig, axes = plt.subplots(len(models), 3, figsize=(15, 3 * len(models)))
    for i, model in enumerate(models):
        plot_mixture(model, data, axes[i, 0], l, u, bins=bins)
        plot_mixture(model, data, axes[i, 1], log=True, l=np.log(l + 0.0001),
                     u=np.log(u), bins=bins)
        plot_probs(model, axes[i, 2], l, u)
    sns.despine(offset=5)
    fig.tight_layout()
    fig.savefig(out_file)


def plot_all_models_bgmm(models, data, l, u, bins, out_file):
    """
    Plot a bunch of BGMMs.

    :param models: list of GMM model objects
    :param data: Ks array
    :param l: lower Ks limit
    :param u: upper Ks limit
    :param bins: number of histogram bins
    :param out_file: output file
    :return: nada
    """
    fig, axes = plt.subplots(len(models), 4, figsize=(20, 3 * len(models)))
    for i, model in enumerate(models):
        plot_mixture(model, data, axes[i, 0], l, u, bins=bins)
        plot_mixture(model, data, axes[i, 1], log=True, l=np.log(l + 0.0001),
                     u=np.log(u), bins=bins)
        plot_probs(model, axes[i, 2], l, u)
        plot_bars_weights(model, axes[i, 3])
    sns.despine(offset=5)
    fig.tight_layout()
    fig.savefig(out_file)


def get_component_probabilities(df, model):
    """
    Get paralogs from a mixture model component.

    :param df: Ks distribution pandas data frame
    :param model: mixture model
    :return: data frame
    """
    df = df.dropna()
    df = df.drop_duplicates(keep='first')
    df['log(Ks)'] = np.log(df['Ks'])
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna()

    p = model.predict_proba(np.array(df['log(Ks)']).reshape(-1, 1))
    order = np.argsort([x[0] for x in model.means_])
    order_dict = {i: order[i] for i in range(len(order))}
    for c in range(len(order)):
        col = 'p_component{}'.format(order_dict[c] + 1)
        df[col] = p[:, c]

    return df
