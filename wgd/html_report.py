#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
import os
import pandas as pd


def write_report_ks(output_dir, species, gene_families, nucleotide_fasta,
                    protein_fasta, counts, mixture='bayesian', kde=True, filename='report.html'):
    """
    Write an `html` report of a Ks distribution analysis.

    :param filename: file name to write report to
    :return: html report written to the specified file
    """
    with open(gene_families, 'r') as f:
        gfs = len(f.readlines())

    with open(nucleotide_fasta, 'r') as f:
        n = len(f.read().split('>'))

    if protein_fasta:
        with open(protein_fasta, 'r') as f:
            p = len(f.read().split('>'))

    else:
        protein_fasta = 'translated from {}'.format(nucleotide_fasta)
        p = n

    Ks = pd.read_csv(os.path.join(output_dir, 'Ks_distribution.csv'), index_col=0)
    Kn = pd.read_csv(os.path.join(output_dir, 'Ka_distribution.csv'), index_col=0)
    w = pd.read_csv(os.path.join(output_dir, 'Omega_distribution.csv'), index_col=0)

    Ks_out = Ks[Ks['Ks'] > 5].shape[0]
    Kn_out = Kn[Kn['Ka'] > 1].shape[0]
    w_out = w[w['Omega'] > 1].shape[0]

    if mixture == 'bayesian':
        mm = BGMM
    elif mixture == 'gaussian':
        mm = GMM
    else:
        mm = 'No mixture modeling performed!'

    if kde:
        kd = KDE
    else:
        kd = "No kernel density estimation performed!"

    with open(os.path.join(output_dir, filename), 'w') as f:
        f.write(HTML_STRING.format(species, gene_families, gfs, nucleotide_fasta, n, protein_fasta, p,
                                   Ks.shape[0], Kn.shape[0], w.shape[0], Ks_out, Kn_out, w_out,
                                   counts, counts, counts, mm, kd))
    return


HTML_STRING = """
<!DOCTYPE html>
<html>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="http://www.w3schools.com/lib/w3.css">
<title>Ks analysis report for {0}</title>

<nav class="w3-sidenav w3-light-green" style="width:15%">
    <div class="w3-container w3-teal w3-green" style="height:75px;">
        <h4><pre> C E D A L I O N </pre></h4>
    </div>
    <a href="#data">Data used</a>
    <a href="#ks_dist"><i>K<sub>S</sub></i> distribution</a>
    <a href="#kn_dist"><i>K<sub>N</sub></i> distribution</a>
    <a href="#w_dist"><i>&omega;</i> distribution</a>
    <a href="#mixture"><i>K<sub>S</sub></i> mixture modeling</a>
    <a href="#kde"><i>K<sub>S</sub></i> Kernel density estimation</a>
</nav>

<div style="margin-left:15%">

    <div class="w3-container w3-teal w3-green" style="height:75px;">
        <h1><i>K<sub>S</sub></i> analysis report</h1>
    </div>

    <div class="w3-container w3-justify" style="width:70%;">

        <p style="font-size:12px;">
            <i><b>Note:</b> For correctly rendering this report you should be connected to the internet.</i>
        </p>
        <a name="data"><h2>Data used</h2></a>
        <h3>Input data</h3>
        <table class="w3-table-all w3-hoverable">
            <tr class="w3-green">
                <td>Data</td>
                <td>File</td>
                <td># entries</td>
            </tr>
            <tr>
                <td>Gene families</td>
                <td><code>{1}</code></td>
                <td>{2}</td>
            </tr>
            <tr>
                <td>Nucleotide sequences</td>
                <td><code>{3}</code></td>
                <td>{4}</td>
            </tr>
            <tr>
                <td>Protein sequences</td>
                <td><code>{5}</code></td>
                <td>{6}</td>
            </tr>
        </table>

        <h3>Distributions</h3>
        <table class="w3-table-all w3-hoverable">
            <tr class="w3-green">
                <td>Metric</td>
                <td># gene families</td>
                <td># gene pairs</td>
                <td>outliers</td>
            </tr>
            <tr>
                <td><i>K<sub>S</sub></i></td>
                <td>{13}</td>
                <td>{7}</td>
                <td>{10}</td>
            </tr>
            <tr>
                <td><i>K<sub>N</sub></i></td>
                <td>{14}</td>
                <td>{8}</td>
                <td>{11}</td>
            </tr>
            <tr>
                <td>&omega;</td>
                <td>{15}</td>
                <td>{9}</td>
                <td>{12}</td>
            </tr>
        </table>
        <p>
            For <i>K<sub>S</sub></i>, outlying pairs are defined as those with <i>K<sub>S</sub></i> > 5. For
            <i>K<sub>N</sub></i> and &omega; the number of pairs with values larger than 1 is reported (in the case of
            &omega; indicative for positive selection).
        </p>

        <a name="ks_dist"><h2><i>K<sub>S</sub></i> distribution</h2></a>
        <p>
            Distribution of the number of synonymous substitutions per synonymous site (<i>K<sub>S</sub></i>) as
            estimated through Maximum likelihood with <code>codeml</code> from the
            <a href="http://abacus.gene.ucl.ac.uk/software/paml.html"><code>PAML</code> package (Yang 2007)</a>.
            The approach for constructing the distribution largey follows Vanneste et al. (2013).
            As <i>K<sub>S</sub></i> saturation may bias the results at high <i>K<sub>S</sub></i> values, the
            <i>K<sub>S</sub></i>  range is restricted from 0.1 to 5.
        </p>
        <img src="Ks_distribution.png" style="width:100%;">

        <a name="kn_dist"><h2><i>K<sub>N</sub></i> distribution</h2></a>
        <p>
            Distribution of the number of nonsynonymous substitutions per nonsynonymous site (<i>K<sub>N</sub></i>)
            as estimated through Maximum likelihood with <code>codeml</code> from the
            <a href="http://abacus.gene.ucl.ac.uk/software/paml.html"><code>PAML</code> package (Yang 2007)</a>.
            Please note that the x-axis scale ranges from 0 to the maximum <i>K<sub>N</sub></i> value smaller than 5 in
            the results.
        </p>
        <img src="Ka_distribution.png" style="width:100%;">

        <a name="w_dist"><h2><i>&omega;</i> distribution</h2></a>
        <p>
            &omega; is the <i>K<sub>N</sub></i>/<i>K<sub>S</sub></i> ratio and is a simplistic measure for selective
            pressure. &omega; = 1 is indicative for neutral evolution (nonsynonymous substitutions are tolerated as
            much as synonymous ones), &omega; < 1 is indicative for purifying selection (nonsynonymous substitutions
            are disfavored) and &omega; > 1 is indicative for positive or diversifying selection (nonsynonymous
            substitutions occur more than synonymous substitutions). Please note that the x-axis scale ranges from 0 to
            the maximum &omega; value in the results.
        </p>
        <img src="Omega_distribution.png" style="width:100%;">

        <a name="mixture"><h2>Mixture modeling</h2></a>
        <p>
            Mixture model for a <i>K<sub>S</sub></i> range of 0.1 to 2. The range is limited at 0.1 to avoid
            incorporation of potential allelic variants and at 2 to avoid influences of <i>K<sub>S</sub></i>
            saturation effects. The vertical dashed lines indicate the geometric mean of each component and can
            serve as a point estimate of a <i>K<sub>S</sub></i> peak of interest. Note that the distribution is
            modeled as a mixture of lognormal components and that a lognormal distributed variable is expected
            to show multiplicative stochastic variation around the geometric mean (Morrison, 2008).
        </p>
        {16}

        <img src="mixture.png" style="width:100%;">

        <a name="kde"><h2>Kernel density estimation</h2></a>
        {17}
        <img src="kde.png" style="width:100%;">
        <br>
        <br>
    </div>
    <div class="w3-container w3-teal" style="text-align:right;">
        <h5><i>Arthur Zwaenepoel (2017)</i> <pre>C E D A L I O N</pre></h5>
    </div>
</div>
</html>
"""


BGMM = """
<p>
    Mixture modeling is performed using the Bayesian Gaussian mixture modeling approach as implemented in
    <code>sklearn.mixture</code>. It has not been thoroughly studied whether mixture modeling is appropriate for
    <i>K<sub>S</sub></i> based age distributions. As Gaussian mixture models tend to overfit for distributions
    that are not clearly mixtures of Gaussians, no good way for model selection in Gaussian mixture modeling
    exists. Usage of information criteria such as AIC (Akaike information criterion) or BIC (Bayesian
    information criterion) tends to result in preference of the model with the most fitted components.
    For these reasons, Mixture models for different numbers of components are shown here, only to aid visual
    interpretation of <i>K<sub>S</sub></i> distributions and to give an idea for point estimates of the peak of
    interest. Bear in mind that mixture models for <i>K<sub>S</sub></i> distributions might be highly
    misleading.
</p>
<p>
    While no explicit model selection is performed, the Bayesian variant of Gaussian mixture models as applied
    here assigns a weight to each component (see the bar plots on the right) which indicates the relative
    importance of each component. This can be in a way thought of as a form of model selection.
    These weights depend on the &gamma; hyperparameter. A larger &gamma; will cause more weight to be assigned
    to fewer components. In these distributions, the influence of the &gamma; parameter was observed to be
    limited, and 0.01 was empirically observed as a good choice.
</p>
"""

GMM = """
<p>
    Mixture modeling is performed using the Gaussian mixture modeling approach as implemented in
    <code>sklearn.mixture</code>. It has not been thoroughly studied whether mixture modeling is appropriate for
    <i>K<sub>S</sub></i> based age distributions. As Gaussian mixture models tend to overfit for distributions
    that are not clearly mixtures of Gaussians, no good way for model selection in Gaussian mixture modeling
    exists. Usage of information criteria such as AIC (Akaike information criterion) or BIC (Bayesian
    information criterion) tends to result in preference of the model with the most fitted components.
    It can be appropriate to look for 'knees' in the AIC curve to determine the most appropriate number
    of components.
</p>
"""

KDE = """
<p>
    Kernel density estimation (KDE) of the <i>K<sub>S</sub></i> distribution for different <i>K<sub>S</sub></i>
    ranges (left from 0.1 to 5, right from 0.1 to 2). KDE is performed with different bandwidths (default 0.01, 0.05,
    0.10, 0.15 and 0.20, which were found empirically to give appropriate results). Selecting the optimal
    bandwidth among these four values is performed using 4-fold cross validation (CV) using maximum likelihood
    as score (on a subset of the data). However, this tends to overfit the data due to a preference for smaller
    bandwidths. Therefore results for different bandwidths are shown, enabling better visual interpretation.
    The plot of the likelihood for the 4-fold CV is shown below, note that the absolute value of the likelihood is
    shown, so lower values correspond to higher likelihood.
</p>
"""