#!/usr/bin/python3.5
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

This hasn't been update in a while and might be broken in some features. I
recommend doing the mixture modeling in a Jupyter notebook using ``sklearn``
instead of using the CLI provide here.
"""
import numpy as np
import peakutils
import logging
import os
import random
from sklearn import mixture
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import scipy.stats as ss
import plumbum as pb
import matplotlib
if not 'DISPLAY' in pb.local.env:
    matplotlib.use('Agg')  # use this backend when no X server
import matplotlib.pyplot as plt
import seaborn as sns


def reflected_kde(data_frame, min_ks, max_ks, out_file):
    """
    Perform Kernel density estimation (KDE) with reflected data.
    """
    df = data_frame.groupby(['Family', 'Node']).mean()
    ks = np.array(df['Ks'])
    ks = ks[ks >= min_ks]
    ks = list(ks[ks <= max_ks])
    ks_reflected = ks + list(-1*np.array(ks))
    ks = np.array(ks)
    
    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.distplot(
        ks_reflected, bins=int(max_ks*50), ax = ax, 
        hist_kws={"rwidth": 0.8, "color": "k", "alpha": 0.2}, 
        color="forestgreen", label="Reflected"
    )
    ax.set_xlim(min_ks, max_ks)
    ax.set_yticks([])
    sns.despine(offset=5, trim=True)
    fig.savefig(out_file, bbox_inches='tight')


def gmm(X, n, **kwargs):
    """
    Compute Gaussian mixtures for different numbers of components
    """
    # fit models with 1 to n components
    N = np.arange(1, n + 1)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        models[i] = mixture.GaussianMixture(
            n_components=N[i], covariance_type='full', max_iter=100, **kwargs
        ).fit(X)

    # compute the AIC and the BIC
    aic = [m.aic(X) for m in models]
    bic = [m.bic(X) for m in models]
    best = models[np.argmin(bic)]
    
    return models, bic, aic, best


def plot_mixture(model, data, ax, l=0.1, u=5, color='black', alpha=0.2, 
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
    return ax


def inspect_aic(aic):
    """
    Evaluate the information criterion results for mixture models.
    """
    im = np.argmin(aic) 
    print("min(AIC) = {:.2f} for model {}".format(aic[im], im+1))
    print("Relative probabilities compared to model {}:".format(im+1))
    print("   /                          \\")
    print("   |      (min(AIC) - AICi)/2 |")
    print("   | p = e                    |")
    print("   \                          /")
    for i, aic_i in enumerate(aic):
        p_i = np.exp((aic[im] - aic_i)/2)
        print(". model{:>4}: p = {:.4f}".format(i+1, p_i))
    print()
        
        
def inspect_bic(bic):
    im = np.argmin(bic)
    l = [
        "0 to  2:   Very weak",
        "2 to  6:    Positive",
        "6 to 10:      Strong",
        "    >10: Very Strong"
    ]
    print("Delta BIC assessment: ")
    print("min(BIC) = {:.2f}".format(bic[im]))
    for i, bic_i in enumerate(bic):
        dbic = bic_i - bic[im]
        j = 0
        if dbic > 2: j = 1
        if dbic > 6: j = 2
        if dbic > 10: j = 3
        print(". model{:>4}: delta(BIC) = {:>8.2f} ({})".format(i+1, dbic, l[j]))
    print()
        
        
def plot_probs(m, ax, l=0.01, u=5, ylab=True):
    """
    Plot posterior component probabilities
    """
    x = np.linspace(l, u, 1000).reshape(-1, 1)
    p = m.predict_proba(np.log(x))
    p = p.cumsum(1).T
    x = np.linspace(l, u, 1000)
    order = tuple(np.argsort(m.means_.reshape((1, -1)))[0])
    alphas = np.linspace(0.2, 1, p.shape[0])[order,]
    for i, array in enumerate(p):
        if i == 0:
            ax.fill_between(
                x, 0, p[i], color='gray', alpha=alphas[i], label=str(i+1))
        elif i == len(alphas) - 1:
            ax.fill_between(
                x, p[i-1], 1, color='gray', alpha=alphas[i], label=str(i+1))
        else:
            ax.fill_between(
                x, p[i-1], p[i], color='gray', alpha=alphas[i], label=str(i+1))
    ax.set_xlim(l, u)
    ax.set_ylim(0, 1)
    ax.set_xlabel('$K_S$')
    if ylab:
        ax.set_ylabel('$P(class|K_S)$')
    ax.legend(frameon=True)
    return ax


def plot_aic(aic, ax):
    ax.plot(np.arange(1,len(aic)+1), aic, color='k', marker='o')
    ax.set_xticks(list(range(1,len(aic)+1)))
    ax.grid(ls=":")
    return ax



# OLD API ======================================================================

def weighted_to_unweighted(data_frame):
    """
    Generate an unweighted histogram from a weighted Ks distribution.

    :param data_frame: Ks distribution data frame
    :return: reweighted Ks distribution
    """
    # drop nan values if present
    logging.debug("Dropping NaN containing rows")
    pre = data_frame.shape[0]
    data_frame = data_frame.dropna()
    na = pre - data_frame.shape[0]
    if na != 0:
        logging.warning("Dropped {} NaN containing rows.".format(na))

    m = min(data_frame['WeightOutliersIncluded'])
    l = []
    for i in range(data_frame.shape[0]):
        times = int(data_frame['WeightOutliersIncluded'].iloc[i] / (m))
        l += [data_frame['Ks'].iloc[i]] * times
    return np.array(l).reshape((-1, 1))


def kernel_density_estimation(
        data_frame, bandwidths=(0.01, 0.05, 0.1, 0.15, 0.2), kernel='gaussian',
        output_dir=None
):
    """
    Kernel density estimation (KDE) for Ks distribution

    :param data_frame: pandas data frame with Ks analysis results
    :param bandwidths: the bandwidths to apply in the KDE
    :param kernel: kernel to aply in the KDE (gaussian, epachenikov, top_hat)
    :param output_dir: output directory
    :return: sklearn.KernelDensity objects
    """
    X = weighted_to_unweighted(data_frame)

    n = len(bandwidths)

    # X = weighted_to_unweighted(data_frame)
    X_full = X[X > 0.1]
    X_full = X_full[X_full <= 5].reshape(-1, 1)
    X_plot_full = np.linspace(0, 5, 1000)[:, np.newaxis]

    X_red = X_full[X_full <= 2].reshape(-1, 1)
    X_plot_red = np.linspace(0, 2, 400)[:, np.newaxis]

    fig = plt.figure(figsize=(20, 6 * (n + 1)))

    for i in range(len(bandwidths)):
        ax = fig.add_subplot(n + 1, 2, 2 * i + 1)

        # full histogram
        ax.set_title("KDE of Ks distribution (Ks < 5), bandwidth = {}".format(
                bandwidths[i]))
        ax.set_xlabel("Binned Ks")
        ax.set_xlim(0.1, 5)
        ax.hist(X_full, bins=100, histtype='stepfilled', color="#82c982",
                alpha=0.4, normed=True)

        kde = KernelDensity(kernel=kernel, bandwidth=bandwidths[i]).fit(X_full)
        log_dens = kde.score_samples(X_plot_full)
        ax.plot(X_plot_full[:, 0], np.exp(log_dens), '-',
                label="kernel = '{0}'".format(kernel), color='#5eba5e')

        indexes = peakutils.indexes(np.exp(log_dens),
                                    thres=0.02 / max(log_dens), min_dist=100)
        peaks = X_plot_full[indexes]

        y_min, y_max = ax.get_ylim()
        ax.vlines(peaks, ymin=y_min, ymax=y_max - 0.1, linestyles='dashed',
                  alpha=0.7)
        for j in range(len(indexes)):
            ax.text(s="{:.2f}".format(peaks[j][0]), x=peaks[j], y=y_max - 0.1)

        ax.legend(loc='upper right')

        # reduced histogram
        ax = fig.add_subplot(n + 1, 2, 2 * i + 2)
        ax.set_title("KDE of Ks distribution (Ks < 2), bandwidth = {}".format(
                bandwidths[i]))
        ax.set_xlabel("Binned Ks")
        ax.set_xlim(0.1, 2)
        ax.hist(X_red, bins=40, histtype='stepfilled', color="#82c982",
                alpha=0.4, normed=True)

        kde = KernelDensity(kernel=kernel, bandwidth=bandwidths[i]).fit(X_red)
        log_dens = kde.score_samples(X_plot_red)
        ax.plot(X_plot_red[:, 0], np.exp(log_dens), '-',
                label="kernel = '{0}'".format(kernel), color='#5eba5e')

        indexes = peakutils.indexes(np.exp(log_dens),
                                    thres=0.02 / max(log_dens), min_dist=100)
        peaks = X_plot_red[indexes]

        y_min, y_max = ax.get_ylim()
        ax.vlines(peaks, ymin=y_min, ymax=y_max - 0.1, linestyles='dashed',
                  alpha=0.7)
        for j in range(len(indexes)):
            ax.text(s="{:.2f}".format(peaks[j][0]), x=peaks[j], y=y_max - 0.1)

        ax.legend(loc='upper right')

    # Choose the 'best' bandwidth based on 4 fold cross-validation on subsample
    X_grid = np.array(random.sample(list(X_full), 2000)).reshape(-1, 1)
    grid = GridSearchCV(KernelDensity(kernel=kernel), {'bandwidth': bandwidths},
                        cv=4)
    grid.fit(np.array(X_grid))

    ax = fig.add_subplot(n + 1, 2, 2 * i + 3)
    ax.set_title('Absolute value of the likelihood for different bandwidths')
    ax.plot(bandwidths, np.abs(grid.cv_results_['mean_test_score']),
            color="#a6d9a6")
    ax.scatter(bandwidths, np.abs(grid.cv_results_['mean_test_score']),
               color="#a6d9a6")
    ax.set_xlabel('Bandwidth')
    ax.set_ylabel('Absolute value of likelihood')

    if output_dir:
        fig.savefig(os.path.join(output_dir, 'kde.png'), bbox_inches='tight')

    else:
        fig.show()


def mixture_model_bgmm(data_frame, n_range=(1, 5), Ks_range=(0.1, 2),
                       gamma=0.01, max_iter=1000,
                       plot=True, log=True, output_dir='./', plot_save=True,
                       output_file='mixture.png',
                       fig_size=None, **kwargs):
    """
    Fit a Bayesian Gaussian mixture model and make some plots.

    :param data_frame: pandas data frame from Ks analysis
    :param n_range: range of numbers of components
    :param Ks_range: Ks value interval to model
    :param gamma: gamma hyperparameter for
        :py:func:`sklearn.mixture.BayesianGaussianMixture`
    :param max_iter: maximum number of iterations
    :param plot: make plots (boolean)
    :param log: fit log-normal components rather than Normal
    :param output_dir: output directory
    :param plot_save: save the plot
    :param output_file: output file name
    :param kwargs: other keyword arguments for
        :py:func:`sklearn.mixture.BayesianGaussianMixture`
    :return: mixture models and plots thereof
    """
    X = weighted_to_unweighted(data_frame)
    X = X[X < Ks_range[1]].reshape((-1, 1))
    X = X[X > Ks_range[0]].reshape((-1, 1))
    if log:
        X = np.log(X)

    models = []
    for n in range(n_range[0], n_range[1] + 1):
        logging.info("Fitting BGMM with {} components".format(n))
        models.append(mixture.BayesianGaussianMixture(
                n_components=n,
                covariance_type='full',
                weight_concentration_prior=gamma,
                max_iter=max_iter,
                weight_concentration_prior_type='dirichlet_distribution',
                **kwargs
        ).fit(X))

    if plot:
        logging.info("Plotting BGMMs")

        if not fig_size:
            fig_size = (18, (n_range[1] - n_range[0] + 1) * 6)
        fig = plt.figure(figsize=fig_size)

        for i in range(len(models)):
            x = np.linspace(Ks_range[0], Ks_range[1], 1000).reshape((-1, 1))
            means = models[i].means_
            varcs = models[i].covariances_
            weights = models[i].weights_

            # plot histogram with fitted components
            ax = fig.add_subplot(n_range[1], 2, 2 * i + 1)
            if log:
                ax.hist(np.exp(X), int(Ks_range[1] * 25), normed=True,
                        rwidth=0.8, color="black", alpha=0.2)
            else:
                ax.hist(X, int(Ks_range[1] * 25), normed=True, rwidth=0.8,
                        color="black", alpha=0.2)

            ax.set_xlim(Ks_range[0], Ks_range[1])

            mix = None
            first = True
            for k in range(len(means)):
                # plot component
                if log:
                    ax.plot(x, ss.lognorm.pdf(
                            x, np.sqrt(varcs[k]),
                            scale=np.exp(means[k])) * weights[k],
                            '--k', color='black', alpha=0.4)
                else:
                    ax.plot(x, ss.norm.pdf(
                            x, loc=means[k],
                            scale=np.sqrt(varcs[k])) * weights[k],
                            '--k', color='black', alpha=0.4)

                # add component to mixture
                if first:
                    if log:
                        mix = ss.lognorm.pdf(
                                x, np.sqrt(varcs[k]),
                                scale=np.exp(means[k])) * weights[k]
                    else:
                        mix = ss.norm.pdf(
                                x, loc=means[k],
                                scale=np.sqrt(varcs[k])) * weights[k]
                    first = False
                else:
                    if log:
                        mix += ss.lognorm.pdf(
                                x, np.sqrt(varcs[k]),
                                scale=np.exp(means[k])) * weights[k]
                    else:
                        mix += ss.norm.pdf(
                                x, loc=means[k],
                                scale=np.sqrt(varcs[k])) * weights[k]

            # plot the mixture
            ax.plot(x, mix, color='black')

            ax.set_xlabel('$K_s$')

            # plot vertical lines for geometric means
            y_min, y_max = ax.get_ylim()

            if log:
                means = np.exp(means)

            ax.vlines(means, ymin=y_min, ymax=y_max, color='black', alpha=0.3)
            for j in range(len(means)):
                ax.text(s="{:.2f}".format(means[j][0]), x=means[j],
                        y=y_max - 0.1)

            ax.set_title(
                    'Fitted components and mixture, $\gamma = {}$'.format(
                            gamma))
            sns.despine(ax=ax, offset=5, trim=True)

            # plot weights and means
            wm = [(weights[j], means[j]) for j in range(len(weights))]
            wm = sorted(wm, key=lambda tup: tup[1])

            ax = fig.add_subplot(n_range[1], 2, 2 * i + 2)
            for m in range(len(wm)):
                ax.bar(m, wm[m][0], color="black", alpha=0.2)
                ax.text(m, wm[m][0], "{0:0.2f}".format(wm[m][1][0]),
                        horizontalalignment='center')

            ax.set_ylabel('Weight')
            ax.set_ylim(0, 1.05)
            ax.set_xticks(list(range(models[i].n_components)))
            ax.set_xticklabels(list(range(1, models[i].n_components + 1)))
            ax.set_xlabel('Component')
            ax.set_title('Weight and mean of each component')
            sns.despine(ax=ax, offset=5, trim=True)

        if plot_save:
            fig.savefig(os.path.join(output_dir, output_file),
                        bbox_inches='tight')

    return models


def mixture_model_gmm(
        data_frame, n=4, Ks_range=(0.1, 2), log=True, plot=True,
        output_dir=None, plot_save=True,
        output_file='mixture.png', fig_size=None,
        **kwargs
):
    """
    Fit mixture models using Gaussian Mixture Modeling for different numbers of
    components and choose best one.

    :param data_frame: pandas data frame from Ks analysis
    :param n: max number of components
    :param Ks_range: Ks range to model
    :param log: fit log-normal components
    :param plot: make plots
    :param output_dir: output directory for plots
    :param plot_save: save plots
    :param output_file: output file name
    :param kwargs: other keyword arguments for ``mixture.GaussianMixture``
    :return: Mixture models and plots
    """
    X = weighted_to_unweighted(data_frame)
    X = X[X < Ks_range[1]].reshape((-1, 1))
    X = X[X > Ks_range[0]].reshape((-1, 1))

    if log:
        X = np.log(X)

    # fit models with 1-n components
    N = np.arange(1, n + 1)
    models = [None for i in range(len(N))]

    for i in range(len(N)):
        logging.info("Fitting GMM with {} components".format(i + 1))
        models[i] = mixture.GaussianMixture(n_components=N[i],
                                            covariance_type='full',
                                            max_iter=100, **kwargs).fit(X)

    # compute the AIC and the BIC
    aic = [m.aic(X) for m in models]
    bic = [m.bic(X) for m in models]

    best = models[np.argmin(bic)]

    logging.info('-' * 24)
    logging.info('{0:<4}{1:>10}{2:>10}'.format('n', 'AIC', 'BIC'))
    logging.info('-' * 24)
    for i in range(n):
        logging.info('{0:<4}{1:^10.2f}{2:>10.2f}'.format(i + 1, aic[i], bic[i]))
    logging.info('-' * 24)

    if plot:
        logging.info("Plotting GMMs")

        if not fig_size:
            fig_size = (20, n * 8)
        fig = plt.figure(figsize=fig_size)

        for i in range(len(models)):
            x = np.linspace(Ks_range[0], Ks_range[1], 1000).reshape((-1, 1))
            means = models[i].means_
            varcs = models[i].covariances_
            weights = models[i].weights_

            # plot histogram with fitted components
            ax = fig.add_subplot(n + 1, 2, 2 * i + 1)
            if log:
                ax.hist(np.exp(X), int(Ks_range[1] * 25), normed=True,
                        color="black", alpha=0.2, rwidth=0.8)
            else:
                ax.hist(X, int(Ks_range[1] * 25), normed=True, color="black",
                        alpha=0.2, rwidth=0.8)
            ax.set_xlim(Ks_range[0], Ks_range[1])

            mix = None
            first = True
            for k in range(len(means)):
                if log:
                    ax.plot(x, ss.lognorm.pdf(
                            x, np.sqrt(varcs[k]),
                            scale=np.exp(means[k])) * weights[k],
                            '--k', color='black', alpha=0.4)
                else:
                    ax.plot(x, ss.norm.pdf(
                            x, loc=means[k],
                            scale=np.sqrt(varcs[k])) * weights[k],
                            '--k', color='black', alpha=0.4)
                if first:
                    if log:
                        mix = ss.lognorm.pdf(
                                x, np.sqrt(varcs[k]),
                                scale=np.exp(means[k])) * weights[k]
                    else:
                        mix = ss.norm.pdf(
                                x, loc=means[k],
                                scale=np.sqrt(varcs[k])) * weights[k]
                    first = False
                else:
                    if log:
                        mix += ss.lognorm.pdf(
                                x, np.sqrt(varcs[k]),
                                scale=np.exp(means[k])) * weights[k]
                    else:
                        mix += ss.norm.pdf(
                                x, loc=means[k],
                                scale=np.sqrt(varcs[k])) * weights[k]

            # plot the mixture
            ax.plot(x, mix, color='black')

            ax.set_xlabel('$K_s$')

            # plot vertical lines for geometric means
            y_min, y_max = ax.get_ylim()

            if log:
                means = np.exp(means)

            ax.vlines(means, ymin=y_min, ymax=y_max, color='black', alpha=0.3)
            for j in range(len(means)):
                ax.text(s="{:.2f}".format(means[j][0]), x=means[j],
                        y=y_max - 0.1)

            sns.despine(ax=ax, offset=5, trim=True)

            # plot weights and means
            wm = [(weights[j], means[j]) for j in range(len(weights))]
            wm = sorted(wm, key=lambda tup: tup[1])

            ax = fig.add_subplot(n + 1, 2, 2 * i + 2)
            for m in range(len(wm)):
                ax.bar(m, wm[m][0], color="black", alpha=0.2)
                ax.text(m, wm[m][0], "{0:0.2f}".format(wm[m][1][0]),
                        horizontalalignment='center')

            ax.set_ylabel('Weight')
            ax.set_ylim(0, 1.05)
            ax.set_xticks(list(range(models[i].n_components)))
            ax.set_xticklabels(list(range(1, models[i].n_components + 1)))
            ax.set_xlabel('Component')
            ax.set_title('Weight and mean of each component')
            sns.despine(ax=ax, offset=5, trim=True)

        ax = fig.add_subplot(n + 1, 2, 2 * i + 3)
        ax.plot(list(range(1, n + 1)), aic, color='black')
        ax.set_xticks(list(range(1, n + 1)))
        ax.set_xlabel("Number of components")
        ax.set_ylabel("AIC")

        ax = fig.add_subplot(n + 1, 2, 2 * i + 4)
        ax.plot(list(range(1, n + 1)), bic, color='black')
        ax.set_xticks(list(range(1, n + 1)))
        ax.set_xlabel("Number of components")
        ax.set_ylabel("BIC")

        if plot_save:
            fig.savefig(os.path.join(output_dir, output_file),
                        bbox_inches='tight')
            plt.close()

    return models, bic, aic, best


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
