#!/usr/bin/python3.5
"""
Arthur Zwaenepoel - 2017
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
        times = int(data_frame['WeightOutliersIncluded'].iloc[i]/(m))
        l += [data_frame['Ks'].iloc[i]]*times
    return np.array(l).reshape((-1,1))


def kernel_density_estimation(data_frame, bandwidths=(0.01, 0.05, 0.1, 0.15, 0.2), kernel='gaussian', output_dir=None):
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
        ax.set_title("KDE of Ks distribution (Ks < 5), bandwidth = {}".format(bandwidths[i]))
        ax.set_xlabel("Binned Ks")
        ax.set_xlim(0.1, 5)
        ax.hist(X_full, bins=100, histtype='stepfilled', color="#82c982", alpha=0.4, normed=True)

        kde = KernelDensity(kernel=kernel, bandwidth=bandwidths[i]).fit(X_full)
        log_dens = kde.score_samples(X_plot_full)
        ax.plot(X_plot_full[:, 0], np.exp(log_dens), '-', label="kernel = '{0}'".format(kernel), color='#5eba5e')

        indexes = peakutils.indexes(np.exp(log_dens), thres=0.02 / max(log_dens), min_dist=100)
        peaks = X_plot_full[indexes]

        y_min, y_max = ax.get_ylim()
        ax.vlines(peaks, ymin=y_min, ymax=y_max-0.1, linestyles='dashed', alpha=0.7)
        for j in range(len(indexes)):
            ax.text(s="{:.2f}".format(peaks[j][0]), x=peaks[j], y=y_max-0.1)

        ax.legend(loc='upper right')

        # reduced histogram
        ax = fig.add_subplot(n + 1, 2, 2 * i + 2)
        ax.set_title("KDE of Ks distribution (Ks < 2), bandwidth = {}".format(bandwidths[i]))
        ax.set_xlabel("Binned Ks")
        ax.set_xlim(0.1, 2)
        ax.hist(X_red, bins=40, histtype='stepfilled', color="#82c982", alpha=0.4, normed=True)

        kde = KernelDensity(kernel=kernel, bandwidth=bandwidths[i]).fit(X_red)
        log_dens = kde.score_samples(X_plot_red)
        ax.plot(X_plot_red[:, 0], np.exp(log_dens), '-', label="kernel = '{0}'".format(kernel), color='#5eba5e')

        indexes = peakutils.indexes(np.exp(log_dens), thres=0.02 / max(log_dens), min_dist=100)
        peaks = X_plot_red[indexes]

        y_min, y_max = ax.get_ylim()
        ax.vlines(peaks, ymin=y_min, ymax=y_max-0.1, linestyles='dashed', alpha=0.7)
        for j in range(len(indexes)):
            ax.text(s="{:.2f}".format(peaks[j][0]), x=peaks[j], y=y_max-0.1)

        ax.legend(loc='upper right')

    # Choose the 'best' bandwidth based on 4 fold cross-validation on subsample
    X_grid = np.array(random.sample(list(X_full), 2000)).reshape(-1, 1)
    grid = GridSearchCV(KernelDensity(kernel=kernel), {'bandwidth': bandwidths}, cv=4)
    grid.fit(np.array(X_grid))

    ax = fig.add_subplot(n + 1, 2, 2 * i + 3)
    ax.set_title('Absolute value of the likelihood for different bandwidths')
    ax.plot(bandwidths, np.abs(grid.cv_results_['mean_test_score']), color="#a6d9a6")
    ax.scatter(bandwidths, np.abs(grid.cv_results_['mean_test_score']), color="#a6d9a6")
    ax.set_xlabel('Bandwidth')
    ax.set_ylabel('Absolute value of likelihood')

    if output_dir:
        fig.savefig(os.path.join(output_dir, 'kde.png'), bbox_inches='tight')

    else:
        fig.show()


def mixture_model_bgmm(data_frame, n_range=(1,5), Ks_range=(0.1, 2), gamma=0.01, max_iter=1000,
                       plot=True, log=True, output_dir='./', plot_save=True, output_file='mixture.png',
                       fig_size=None, **kwargs):
    """
    Fit a Bayesian Gaussian mixture model and make some plots.

    :param data_frame: pandas data frame from Ks analysis
    :param n_range: range of numbers of components
    :param Ks_range: Ks value interval to model
    :param gamma: gamma hyperparameter for sklearn.mixture.BayesianGaussianMixture
    :param max_iter: maximum number of iterations
    :param plot: make plots (boolean)
    :param log: fit log-normal components rather than Normal
    :param output_dir: output directory
    :param plot_save: save the plot
    :param output_file: output file name
    :param kwargs: other keyword arguments for :py:func:`sklearn.mixture.BayesianGaussianMixture`
    :return: mixture models and plots thereof
    """
    X = weighted_to_unweighted(data_frame)
    X = X[X < Ks_range[1]].reshape((-1, 1))
    X = X[X > Ks_range[0]].reshape((-1, 1))
    if log:
        X = np.log(X)

    models = []
    for n in range(n_range[0], n_range[1]+1):
        logging.info("Fitting BGMM with {} components".format(n))
        models.append(mixture.BayesianGaussianMixture(n_components=n, covariance_type='full',
                                                      weight_concentration_prior=gamma, max_iter=max_iter,
                                                      weight_concentration_prior_type='dirichlet_distribution',
                                                      **kwargs).fit(X))

    if plot:
        logging.info("Plotting BGMMs")

        if not fig_size:
            fig_size = (20, (n_range[1]-n_range[0]+1) * 6)
        fig = plt.figure(figsize=fig_size)

        for i in range(len(models)):
            x = np.linspace(Ks_range[0], Ks_range[1], 1000).reshape((-1, 1))
            means = models[i].means_
            varcs = models[i].covariances_
            weights = models[i].weights_

            # plot histogram with fitted components
            ax = fig.add_subplot(n_range[1], 2, 2 * i + 1)
            if log:
                ax.hist(np.exp(X), int(Ks_range[1] * 25), normed=True, rwidth=0.8, color="black", alpha=0.2)
            else:
                ax.hist(X, int(Ks_range[1] * 25), normed=True, rwidth=0.8, color="black", alpha=0.2)

            ax.set_xlim(Ks_range[0], Ks_range[1])

            mix = None
            first = True
            for k in range(len(means)):
                # plot component
                if log:
                    ax.plot(x, ss.lognorm.pdf(x, np.sqrt(varcs[k]), scale=np.exp(means[k])) * weights[k],
                            '--k', color='black', alpha=0.4)
                else:
                    ax.plot(x, ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k],
                            '--k', color='black', alpha=0.4)

                # add component to mixture
                if first:
                    if log:
                        mix = ss.lognorm.pdf(x, np.sqrt(varcs[k]), scale=np.exp(means[k])) * weights[k]
                    else:
                        mix = ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k]
                    first = False
                else:
                    if log:
                        mix += ss.lognorm.pdf(x, np.sqrt(varcs[k]), scale=np.exp(means[k])) * weights[k]
                    else:
                        mix += ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k]

            # plot the mixture
            ax.plot(x, mix, color='black')

            ax.set_xlabel('$K_s$')

            # plot vertical lines for geometric means
            y_min, y_max = ax.get_ylim()

            if log:
                means = np.exp(means)

            ax.vlines(means, ymin=y_min, ymax=y_max, color='black', alpha=0.3)
            for j in range(len(means)):
                ax.text(s="{:.2f}".format(means[j][0]), x=means[j], y=y_max - 0.1)

            ax.set_title('Fitted components and mixture, $\gamma = {}$'.format(gamma))
            sns.despine(ax=ax, offset=5, trim=True)

            # plot weights and means
            wm = [(weights[j], means[j]) for j in range(len(weights))]
            wm = sorted(wm, key=lambda tup: tup[1])

            ax = fig.add_subplot(n_range[1], 2, 2 * i + 2)
            for m in range(len(wm)):
                ax.bar(m, wm[m][0], color="black", alpha=0.2)
                ax.text(m, wm[m][0], "{0:0.2f}".format(wm[m][1][0]), horizontalalignment='center')

            ax.set_ylabel('Weight')
            ax.set_ylim(0, 1.05)
            ax.set_xticks(list(range(models[i].n_components)))
            ax.set_xticklabels(list(range(1, models[i].n_components + 1)))
            ax.set_xlabel('Component')
            ax.set_title('Weight and mean of each component')
            sns.despine(ax=ax, offset=5, trim=True)

        if plot_save:
            fig.savefig(os.path.join(output_dir, output_file), bbox_inches='tight')

    return models


def mixture_model_gmm(data_frame, n=4, Ks_range=(0.1, 2), log=True, plot=True,
                      output_dir=None, plot_save=True, output_file='mixture.png', fig_size=None,
                      **kwargs):
    """
    Fit mixture models using Gaussian Mixture Modeling for different numbers of components and choose best one.

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
        logging.info("Fitting GMM with {} components".format(i))
        models[i] = mixture.GaussianMixture(n_components=N[i], covariance_type='full',
                                            max_iter=100, **kwargs).fit(X)

    # compute the AIC and the BIC
    aic = [m.aic(X) for m in models]
    bic = [m.bic(X) for m in models]

    best = models[np.argmin(bic)]

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
                ax.hist(np.exp(X), int(Ks_range[1] * 25), normed=True, color="black", alpha=0.2, rwidth=0.8)
            else:
                ax.hist(X, int(Ks_range[1] * 25), normed=True, color="black", alpha=0.2, rwidth=0.8)
            ax.set_xlim(Ks_range[0], Ks_range[1])

            mix = None
            first = True
            for k in range(len(means)):
                if log:
                    ax.plot(x, ss.lognorm.pdf(x, np.sqrt(varcs[k]), scale=np.exp(means[k])) * weights[k],
                            '--k', color='black', alpha=0.4)
                else:
                    ax.plot(x, ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k],
                            '--k', color='black', alpha=0.4)
                if first:
                    if log:
                        mix = ss.lognorm.pdf(x, np.sqrt(varcs[k]), scale=np.exp(means[k])) * weights[k]
                    else:
                        mix = ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k]
                    first = False
                else:
                    if log:
                        mix += ss.lognorm.pdf(x, np.sqrt(varcs[k]), scale=np.exp(means[k])) * weights[k]
                    else:
                        mix += ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k]

            # plot the mixture
            ax.plot(x, mix, color='black')

            ax.set_xlabel('$K_s$')

            # plot vertical lines for geometric means
            y_min, y_max = ax.get_ylim()

            if log:
                means = np.exp(means)

            ax.vlines(means, ymin=y_min, ymax=y_max, color='black', alpha=0.3)
            for j in range(len(means)):
                ax.text(s="{:.2f}".format(means[j][0]), x=means[j], y=y_max - 0.1)

            sns.despine(ax=ax, offset=5, trim=True)

            # plot weights and means
            wm = [(weights[j], means[j]) for j in range(len(weights))]
            wm = sorted(wm, key=lambda tup: tup[1])

            ax = fig.add_subplot(n + 1, 2, 2 * i + 2)
            for m in range(len(wm)):
                ax.bar(m, wm[m][0], color="black", alpha=0.2)
                ax.text(m, wm[m][0], "{0:0.2f}".format(wm[m][1][0]), horizontalalignment='center')

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
            fig.savefig(os.path.join(output_dir, output_file), bbox_inches='tight')

    return models, bic, aic, best