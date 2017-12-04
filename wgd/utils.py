#!/usr/bin/python3.5
"""
Arthur Zwaenepoel - 2017
"""
import os
import logging
import re
import warnings
import random
import numpy as np
import json
import subprocess
import uuid
from progressbar import ProgressBar
from numpy import mean, std


def can_i_run_software(software):
    """
    Test if external software is executable

    :param software: list or string of executable(s)
    :return: 1 (failure) or 0 (success)
    """
    if type(software) == str:
        software = [software]
    ex = 0
    for s in software:
        if s == 'codeml':
            tmp_file = str(uuid.uuid4())
            with open(tmp_file, 'w') as o:
                o.write('test')
            command = ['codeml', tmp_file]
        else:
            command = [s, '-h']
        try:
            subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            logging.error('{} executable not found!'.format(s))
            ex = 1
        if s == 'codeml':
            os.remove(tmp_file)
    return ex


def get_gfs_for_species(gene_family_dict, gene_pattern):
    """
    Get non-singleton gene families for a species of interest
    :param gene_family_dict: gene family dictionary
    :param species: species of interest
    :return: dictionairy with gene families
    """
    gene_pattern = re.compile(gene_pattern)
    gf_dict = {}

    for key in gene_family_dict.keys():
        if gene_pattern.search(' '.join(gene_family_dict[key])) is not None:
            gf_dict[key] = gene_family_dict[key]

    return gf_dict


def get_sequences(paralog_dict, sequences):
    """
    Fetch sequences from a fasta file or sequence dict and put them in a two level dictionairy
    {gene_family: {gene: seq, gene: seq, ...}, ...}
    :param paralog_dict:
    :return: two-level dictionairy
    """

    if not type(sequences) == dict:
        sequences = read_fasta(sequences, split_on_pipe=True)

    paralog_sequence_dict = {}

    for family in paralog_dict:
        paralog_sequence_dict[family] = {}
        for gene in paralog_dict[family]:
            if gene not in sequences.keys():
                warnings.warn("Gene {} in gene families but not in protein fasta!".format(gene))
            else:
                paralog_sequence_dict[family][gene] = sequences[gene]

    return paralog_sequence_dict


def process_gene_families(gene_family_file, ignore_prefix=False):
    """
    Processes a raw gene family file as e.g. from OrthoMCL into a generic dictionary structure
    OrthoMCL raw file consists of one gene family per line, including tab separated gene IDs,
    without gene family ID.
    """
    gene_family_dict = {}
    ID = 1

    with open(gene_family_file, 'r') as f:
        for line in f:
            genes = line.strip().split("\t")
            if ignore_prefix:
                if '|' in genes[0]:
                    genes = [gene.split('|')[1] for gene in genes]
            gene_family_dict["GF_{:06d}".format(ID)] = genes
            ID += 1

    return gene_family_dict


def check_dirs(tmp_dir, output_dir, prompt, preserve):
    """
    Check directories needed
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param prompt: prompt for overwrites (boolean)?
    :param preserve: preserve MSA files (boolean)?
    :return: nothing
    """
    # Check the tmp directory
    if tmp_dir:
        if os.path.exists(tmp_dir):
            if prompt:
                overwrite = input("tmp directory {} already exists. Overwrite? [y/n]: ".format(tmp_dir))
                if overwrite == 'y':
                    os.system('rm -r {}'.format(os.path.join(tmp_dir)))
                else:
                    print('EXIT')
                    return
            else:
                os.system('rm -r {}'.format(os.path.join(tmp_dir)))
        os.mkdir(tmp_dir)

    # Check the output directory
    if output_dir:
        if os.path.exists(output_dir):
            if prompt:
                overwrite = input("Output directory {} already exists. Overwrite? [y/n]: ".format(output_dir))
                if overwrite == 'y':
                    os.system('rm -r {}'.format(os.path.join(output_dir)))
                else:
                    print('EXIT')
                    return
            else:
                os.system('rm -r {}'.format(os.path.join(tmp_dir)))
        os.mkdir(output_dir)

    # If data should be preserved, make/clean directories
    if preserve:
        if 'msa' not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir, 'msa'))

        if 'codeml' not in os.listdir(output_dir):
            os.mkdir(os.path.join(output_dir, 'codeml'))


def read_fasta(fasta_file, prefix=None, split_on_pipe=True, split_on_whitespace=True, raw=False):
    """
    Generic fastafile reader
    Returns a dictionairy {ID: sequence, ID: sequence, ...}.
    """
    sequence_dict = {}

    with open(fasta_file, 'r') as f:
        if raw:
            return f.read()

        genes = f.read().split(">")
        for gene in genes:
            ID = gene.split("\n")[0]
            if ID != '':
                if split_on_pipe:
                    ID = ID.split("|")[0].strip()
                if split_on_whitespace:
                    ID = ID.split()[0]
                if prefix and prefix != '':
                    ID = prefix + '|' + ID
                sequence = "".join(gene.split("\n")[1:])
                sequence = sequence.replace('*', '')
                sequence_dict[ID] = sequence

    if '' in sequence_dict.keys():
        del sequence_dict['']
    return sequence_dict


def get_paralogs_fasta(input_fasta, selected_paralogs, output_fasta):
    """
    Get the fasta file associated with the paralogs in a slice from a Ks distribution data frame

    :param input_fasta: fasta file
    :param selected_paralogs: Data frame (slice)
    :param output_fasta: output fasta file
    :return: nada
    """
    seqs = read_fasta(input_fasta)
    genes = set(selected_paralogs['Paralog1']) | set(selected_paralogs['Paralog2'])
    with open(output_fasta, 'w') as f:
        for gene in list(genes):
            ks_values = list(selected_paralogs[selected_paralogs['Paralog1'] == gene]['Ks'])
            ks_values += list(selected_paralogs[selected_paralogs['Paralog2'] == gene]['Ks'])
            ks_value, var = mean(ks_values), std(ks_values)
            if gene in seqs.keys():
                f.write('>{0} mean(Ks)={1:.5f};std(Ks)={2:.5f}\n{3}\n'.format(
                    gene, float(ks_value), float(var), seqs[gene]))
            else:
                logging.warning('gene {} not found in fasta file!'.format(gene))


def translate_cds(sequence_dict):
    """
    Just another CDS to protein translater

    :param sequence_dict: dictionary with gene IDs and CDS sequences
    :return: dictionary with gene IDs and proteins sequences
    """
    aa_dict = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
    }
    protein_dict = {}

    with ProgressBar(max_value=len(sequence_dict.keys())+1) as pb:
        j = 0
        for key, val in sequence_dict.items():
            j += 1
            aa_seq = ''
            for i in range(0,len(val),3):
                if val[i:i+3] not in aa_dict.keys():
                    logging.debug('Invalid codon {0:>3} in sequence {1}'.format(val[i:i+3], key))
                else:
                    aa_seq += aa_dict[val[i:i+3]]
            protein_dict[key] = aa_seq
            pb.update(j)

    return protein_dict


def write_fasta(seq_dict, output_file):
    """
    Write a sequence dictionary to a fasta file

    :param seq_dict:
    :param output_file:
    """
    with open(output_file, 'w') as o:
        for key, val in seq_dict.items():
            o.write('>' + key + '\n')
            o.write(val + '\n')


def filter_one_vs_one_families(gene_families, s1, s2):
    """
    Filter one-vs-one ortholog containing families for two given species.

    :param gene_families:
    :param s1:
    :param s2:
    :return:
    """
    to_delete = []
    for key, val in gene_families.items():
        count = 0
        for gene in val:
            prefix=gene.split('|')[0]
            if prefix == s1 or prefix == s2:
                count += 1
        if count != 2:
            to_delete.append(key)
    for k in to_delete:
        del gene_families[k]
    return gene_families


def get_number_of_sp(genes):
    """
    Get the number of unique species in a list of gene IDs.
    Will be approximate since it is based on the assumption that
    the leading non-digit part identifies the species
    (which is not always the case e.g. mitochondrial genes etc.)
    """
    p = re.compile('\D*')
    gene_set = set()
    for gene in genes:
        m = p.match(gene)
        if m:
            gene_set.add(m.group())
    return len(list(gene_set))


def check_genes(genes, ids):
    for gene in genes:
        for id in ids:
            if gene.startswith(id):
                return True


def _random_color():
    """
    Generate a random hex color
    """
    def r(): return random.randint(0, 255)
    return '#%02X%02X%02X' % (r(), r(), r())


class Genome:
    """
    Class that represents a structural annotation.
    Collects several nice data structures for a genome and parsers for various
    genomic data file formats (e.g. gff, fasta, ...)
    """

    def __init__(self):
        """
        Genome.genome: dictionary with a full representation
        Genome.gene_lists: dictionary with ordered lists of genes per
        chromosome (as for I-ADHoRe)
        """
        self.parent_file = None
        self.genome = {}
        self.gene_lists = {}
        self.colors = {}

    def parse_plaza_gff(self, gff_file, keyword='mRNA', id_string='Parent'):
        """
        Parse a PLAZA annotation file into a genome dictionary

        :param gff_file: input gff (PLAZA style)
        :param keyword: keyword for elements to parse out
        """
        self.parent_file = gff_file

        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')

                if line[2] == keyword:
                    chromosome = line[0]
                    start = line[3]
                    stop = line[4]
                    orientation = line[6]
                    gene_l = line[8].split(';')
                    gene_dict = {x.split('=')[0]: x.split('=')[1] for x in gene_l if len(x.split('=')) == 2}

                    if chromosome not in self.genome:
                        self.genome[chromosome] = {}
                        self.gene_lists[chromosome] = []
                        self.colors[chromosome] = _random_color()

                    self.genome[chromosome][gene_dict[id_string]] = {
                        'orientation': orientation, 'start': start, 'stop': stop}
                    self.gene_lists[chromosome].append(
                        (gene_dict[id_string], orientation, start, stop))
        return

    def karyotype_json(self, out_file='genome.json'):
        """
        Generate karyotype data file in json format (as per Circos.js/d3.js)
        """
        karyotype = []
        for chrom in self.gene_lists.keys():
            # approximate chromosome length
            coordinates = [int(x[2]) for x in self.gene_lists[chrom]]
            coordinates += [int(x[3]) for x in self.gene_lists[chrom]]
            length = max(coordinates) - min(coordinates)
            karyotype.append({'id': chrom, 'label': chrom,
                              'color': self.colors[chrom], 'len': length})

        if out_file:
            with open(out_file, 'w') as f:
                json.dump(karyotype, f)

        else:
            return json.dumps(karyotype)


from scipy.spatial.distance import cdist


class gaussian_kde(object):
    """
    from: https://stackoverflow.com/questions/27623919/weighted-gaussian-kernel-density-estimation-in-python

    Representation of a kernel-density estimate using Gaussian kernels.

    Kernel density estimation is a way to estimate the probability density
    function (PDF) of a random variable in a non-parametric way.
    `gaussian_kde` works for both uni-variate and multi-variate data.   It
    includes automatic bandwidth determination.  The estimation works best for
    a unimodal distribution; bimodal or multi-modal distributions tend to be
    oversmoothed.

    Parameters
    ----------
    dataset : array_like
        Datapoints to estimate from. In case of univariate data this is a 1-D
        array, otherwise a 2-D array with shape (# of dims, # of data).
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a scalar,
        this will be used directly as `kde.factor`.  If a callable, it should
        take a `gaussian_kde` instance as only parameter and return a scalar.
        If None (default), 'scott' is used.  See Notes for more details.
    weights : array_like, shape (n, ), optional, default: None
        An array of weights, of the same shape as `x`.  Each value in `x`
        only contributes its associated weight towards the bin count
        (instead of 1).

    Attributes
    ----------
    dataset : ndarray
        The dataset with which `gaussian_kde` was initialized.
    d : int
        Number of dimensions.
    n : int
        Number of datapoints.
    neff : float
        Effective sample size using Kish's approximation.
    factor : float
        The bandwidth factor, obtained from `kde.covariance_factor`, with which
        the covariance matrix is multiplied.
    covariance : ndarray
        The covariance matrix of `dataset`, scaled by the calculated bandwidth
        (`kde.factor`).
    inv_cov : ndarray
        The inverse of `covariance`.

    Methods
    -------
    kde.evaluate(points) : ndarray
        Evaluate the estimated pdf on a provided set of points.
    kde(points) : ndarray
        Same as kde.evaluate(points)
    kde.pdf(points) : ndarray
        Alias for ``kde.evaluate(points)``.
    kde.set_bandwidth(bw_method='scott') : None
        Computes the bandwidth, i.e. the coefficient that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        .. versionadded:: 0.11.0
    kde.covariance_factor : float
        Computes the coefficient (`kde.factor`) that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        The default is `scotts_factor`.  A subclass can overwrite this method
        to provide a different method, or set it through a call to
        `kde.set_bandwidth`.

    Notes
    -----
    Bandwidth selection strongly influences the estimate obtained from the KDE
    (much more so than the actual shape of the kernel).  Bandwidth selection
    can be done by a "rule of thumb", by cross-validation, by "plug-in
    methods" or by other means; see [3]_, [4]_ for reviews.  `gaussian_kde`
    uses a rule of thumb, the default is Scott's Rule.

    Scott's Rule [1]_, implemented as `scotts_factor`, is::

        n**(-1./(d+4)),

    with ``n`` the number of data points and ``d`` the number of dimensions.
    Silverman's Rule [2]_, implemented as `silverman_factor`, is::

        (n * (d + 2) / 4.)**(-1. / (d + 4)).

    Good general descriptions of kernel density estimation can be found in [1]_
    and [2]_, the mathematics for this multi-dimensional implementation can be
    found in [1]_.

    References
    ----------
    .. [1] D.W. Scott, "Multivariate Density Estimation: Theory, Practice, and
           Visualization", John Wiley & Sons, New York, Chicester, 1992.
    .. [2] B.W. Silverman, "Density Estimation for Statistics and Data
           Analysis", Vol. 26, Monographs on Statistics and Applied Probability,
           Chapman and Hall, London, 1986.
    .. [3] B.A. Turlach, "Bandwidth Selection in Kernel Density Estimation: A
           Review", CORE and Institut de Statistique, Vol. 19, pp. 1-33, 1993.
    .. [4] D.M. Bashtannyk and R.J. Hyndman, "Bandwidth selection for kernel
           conditional density estimation", Computational Statistics & Data
           Analysis, Vol. 36, pp. 279-298, 2001.

    Examples
    --------
    Generate some random two-dimensional data:

    >>> from scipy import stats
    >>> def measure(n):
    >>>     "Measurement model, return two coupled measurements."
    >>>     m1 = np.random.normal(size=n)
    >>>     m2 = np.random.normal(scale=0.5, size=n)
    >>>     return m1+m2, m1-m2

    >>> m1, m2 = measure(2000)
    >>> xmin = m1.min()
    >>> xmax = m1.max()
    >>> ymin = m2.min()
    >>> ymax = m2.max()

    Perform a kernel density estimate on the data:

    >>> X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    >>> positions = np.vstack([X.ravel(), Y.ravel()])
    >>> values = np.vstack([m1, m2])
    >>> kernel = stats.gaussian_kde(values)
    >>> Z = np.reshape(kernel(positions).T, X.shape)

    Plot the results:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
    ...           extent=[xmin, xmax, ymin, ymax])
    >>> ax.plot(m1, m2, 'k.', markersize=2)
    >>> ax.set_xlim([xmin, xmax])
    >>> ax.set_ylim([ymin, ymax])
    >>> plt.show()

    """

    def __init__(self, dataset, bw_method=None, weights=None):
        self.dataset = np.atleast_2d(dataset)
        if not self.dataset.size > 1:
            raise ValueError("`dataset` input should have multiple elements.")
        self.d, self.n = self.dataset.shape

        if weights is not None:
            self.weights = weights / np.sum(weights)
        else:
            self.weights = np.ones(self.n) / self.n

        # Compute the effective sample size
        # http://surveyanalysis.org/wiki/Design_Effects_and_Effective_Sample_Size#Kish.27s_approximate_formula_for_computing_effective_sample_size
        self.neff = 1.0 / np.sum(self.weights ** 2)

        self.set_bandwidth(bw_method=bw_method)

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError : if the dimensionality of the input points is different than
                     the dimensionality of the KDE.

        """
        points = np.atleast_2d(points)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = np.reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                                                                              self.d)
                raise ValueError(msg)

        # compute the normalised residuals
        chi2 = cdist(points.T, self.dataset.T, 'mahalanobis', VI=self.inv_cov) ** 2
        # compute the pdf
        result = np.sum(np.exp(-.5 * chi2) * self.weights, axis=1) / self._norm_factor

        return result

    __call__ = evaluate

    def scotts_factor(self):
        return np.power(self.neff, -1. / (self.d + 4))

    def silverman_factor(self):
        return np.power(self.neff * (self.d + 2.0) / 4.0, -1. / (self.d + 4))

    #  Default method to calculate bandwidth, can be overwritten by subclass
    covariance_factor = scotts_factor

    def set_bandwidth(self, bw_method=None):
        """Compute the estimator bandwidth with given method.

        The new bandwidth calculated after a call to `set_bandwidth` is used
        for subsequent evaluations of the estimated density.

        Parameters
        ----------
        bw_method : str, scalar or callable, optional
            The method used to calculate the estimator bandwidth.  This can be
            'scott', 'silverman', a scalar constant or a callable.  If a
            scalar, this will be used directly as `kde.factor`.  If a callable,
            it should take a `gaussian_kde` instance as only parameter and
            return a scalar.  If None (default), nothing happens; the current
            `kde.covariance_factor` method is kept.

        Notes
        -----
        .. versionadded:: 0.11

        Examples
        --------
        >>> x1 = np.array([-7, -5, 1, 4, 5.])
        >>> kde = stats.gaussian_kde(x1)
        >>> xs = np.linspace(-10, 10, num=50)
        >>> y1 = kde(xs)
        >>> kde.set_bandwidth(bw_method='silverman')
        >>> y2 = kde(xs)
        >>> kde.set_bandwidth(bw_method=kde.factor / 3.)
        >>> y3 = kde(xs)

        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.plot(x1, np.ones(x1.shape) / (4. * x1.size), 'bo',
        ...         label='Data points (rescaled)')
        >>> ax.plot(xs, y1, label='Scott (default)')
        >>> ax.plot(xs, y2, label='Silverman')
        >>> ax.plot(xs, y3, label='Const (1/3 * Silverman)')
        >>> ax.legend()
        >>> plt.show()

        """
        if bw_method is None:
            pass
        elif bw_method == 'scott':
            self.covariance_factor = self.scotts_factor
        elif bw_method == 'silverman':
            self.covariance_factor = self.silverman_factor
        elif np.isscalar(bw_method) and not isinstance(bw_method, string_types):
            self._bw_method = 'use constant'
            self.covariance_factor = lambda: bw_method
        elif callable(bw_method):
            self._bw_method = bw_method
            self.covariance_factor = lambda: self._bw_method(self)
        else:
            msg = "`bw_method` should be 'scott', 'silverman', a scalar " \
                  "or a callable."
            raise ValueError(msg)

        self._compute_covariance()

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor().
        """
        self.factor = self.covariance_factor()
        # Cache covariance and inverse covariance of the data
        if not hasattr(self, '_data_inv_cov'):
            # Compute the mean and residuals
            _mean = np.sum(self.weights * self.dataset, axis=1)
            _residual = (self.dataset - _mean[:, None])
            # Compute the biased covariance
            self._data_covariance = np.atleast_2d(np.dot(_residual * self.weights, _residual.T))
            # Correct for bias (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_covariance)
            self._data_covariance /= (1 - np.sum(self.weights ** 2))
            self._data_inv_cov = np.linalg.inv(self._data_covariance)

        self.covariance = self._data_covariance * self.factor ** 2
        self.inv_cov = self._data_inv_cov / self.factor ** 2
        self._norm_factor = np.sqrt(np.linalg.det(2 * np.pi * self.covariance))  # * self.n
