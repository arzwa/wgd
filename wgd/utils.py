# TODO: clean up, currently obsolete, but we may wish to recycle stuff...
import os
import logging
import re
import warnings
import random
import numpy as np
import json
import subprocess
from progressbar import ProgressBar
from numpy import mean, std
from scipy.spatial.distance import cdist


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
        # codeml needs input otherwise it prompts the user for input, so a dummy
        # file is created
        if s == 'codeml':
            tmp_file = uniq_id()
            with open(tmp_file, 'w') as o:
                o.write('test')
            command = ['codeml', tmp_file]
        elif s == 'prank':
            command = [s, '--help']
        elif s in ['blastp', 'makeblastdb', 'blast', 'muscle', 'i-adhore']:
            command = [s, '-version']
        else:
            command = [s, '--version']
        try:
            sp = subprocess.run(command, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            logging.info(sp.stdout.decode("utf-8").strip())
        except FileNotFoundError:
            logging.error('{} executable not found!'.format(s))
            ex = 1

        # remove the dummy file
        if s == 'codeml':
            logging.info('codeml found')
            os.remove(tmp_file)
            subprocess.run(['rm', 'rub', 'rst1', 'rst'])
    return ex


def log_subprocess(program, process):
    """
    Log output from a subprocess call to debug log stream

    :param program: program name
    :param process: completed subprocess object
    """
    logging.debug('{} stdout:\n'.format(program) +
                  process.stdout.decode('utf-8'))
    logging.debug('{} stderr:\n'.format(program) +
                  process.stderr.decode('utf-8'))


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
    Fetch sequences from a fasta file or sequence dict and put them in a two
    level dictionairy {gene_family: {gene: seq, gene: seq, ...}, ...}

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
                warnings.warn("Gene {} in gene families but not in protein "
                              "fasta!".format(gene))
            else:
                paralog_sequence_dict[family][gene] = sequences[gene]

    return paralog_sequence_dict


def process_gene_families(gene_family_file, ignore_prefix=False):
    """
    Processes a raw gene family file as e.g. from MCL output into a generic
    dictionary structure. MCL raw file consists of one gene family per line,
    including tab separated gene IDs, (without gene family ID !).

    Example::

        gene1   gene2   gene3   gene4   gene5
        gene6   gene7   gene8
        gene9   gene10  gene11
        gene12

    :param gene_family_file: file in the right raw gene family format
        (see above)
    :param ignore_prefix: ignore prefixes (boolean), if the gene contains a '|'
        character in default mode the part preceding the '|' is trimmed of,
        assuming the second part is the gene ID. If ``ignore prefix = True``
        this behavior is suppressed. So with ``ignore_prefix = False``
        ``ath|AT1G10000`` becomes ``AT1G10000`` after processing, with
        ``ignore_prefix = True`` it remains ``ath|AT1G10000``.
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


def uniq_id():
    """
    Get a unique ID which is not crazy long

    :return: string
    """
    from time import time
    return str(hex(int(time()*10000000))[2:])


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
                overwrite = input("tmp directory {} already exists. Overwrite? "
                                  "[y/n]: ".format(tmp_dir))
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
                overwrite = input("Output directory {} already exists. "
                                  "Overwrite? [y/n]: ".format(output_dir))
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


def read_fasta(fasta_file, prefix=None, split_on_pipe=False,
               split_on_whitespace=True, raw=False):
    """
    Generic fasta file reader. Returns a dictionairy ``{ID: sequence, ID:
    sequence, ...}``.

    :param fasta_file: fasta file
    :param prefix: prefix to add to gene IDs
    :param split_on_pipe: boolean, split gene IDs on '|' character
    :param split_on_whitespace: boolean, split gene IDs on whitespace
    :param raw: boolean, return raw fasta file (why would you want this?)
    :return: sequence dictionary
    """
    invalid_chars = [":", ";", "?"]
    sequence_dict = {}

    with open(fasta_file, 'r') as f:
        if raw:
            return f.read()

        genes = f.read().split(">")
        for gene in genes:
            ID = gene.split("\n")[0]
            for char in invalid_chars:
                if char in ID:
                    raise ValueError(
                        "Gene ID '{0}' contains illegal character '{1}'".format(ID, char))
            if ID != '':
                if split_on_pipe:
                    ID = ID.split("|")[0].strip()
                if split_on_whitespace:
                    ID = ID.split()[0]
                if prefix and prefix != '':
                    ID = prefix + '|' + ID
                sequence = "".join(gene.split("\n")[1:])
                sequence = sequence.replace('*', '')
                sequence_dict[ID] = sequence.upper()

    if '' in sequence_dict.keys():
        del sequence_dict['']
    return sequence_dict


def get_paralogs_fasta(input_fasta, selected_paralogs, output_fasta,
                       pairs=False):
    """
    Get the fasta file associated with the paralogs in a slice from a Ks
    distribution data frame

    :param input_fasta: fasta file
    :param selected_paralogs: Data frame (slice)
    :param output_fasta: output fasta file
    :return: nada
    """
    seqs = read_fasta(input_fasta)

    if not pairs:
        genes = set(selected_paralogs['Paralog1']) | \
                set(selected_paralogs['Paralog2'])
        with open(output_fasta, 'w') as f:
            for gene in list(genes):
                ks_values = list(
                        selected_paralogs[
                            selected_paralogs['Paralog1'] == gene]['Ks']
                )
                ks_values += list(
                        selected_paralogs[
                            selected_paralogs['Paralog2'] == gene]['Ks']
                )
                ks_value, var = mean(ks_values), std(ks_values)
                if gene in seqs.keys():
                    f.write('>{0} mean(Ks)={1:.5f};std(Ks)={2:.5f}\n{3}\n'
                            ''.format(
                            gene, float(ks_value), float(var), seqs[gene]))
                else:
                    logging.warning(
                            'Gene {} not found in fasta file!'.format(gene))

    if pairs:
        for row in selected_paralogs.index:
            p1, p2 = selected_paralogs.loc[row]['Paralog1'], \
                     selected_paralogs.loc[row]['Paralog2']
            if p1 in seqs.keys() and p2 in seqs.keys():
                with open('{0}_{1}_{2:.3f}_{3}'.format(
                        p1, p2, selected_paralogs.loc[row]['Ks'], output_fasta),
                        'w') as o:
                    o.write('>{0}\n{1}\n>{2}\n{3}'.format(
                            p1, seqs[p1], p2, seqs[p2]))
            else:
                logging.warning('Gene not found in fasta file!')

    return


def translate_cds(sequence_dict, skip_invalid=False):
    """
    Just another CDS to protein translater. Will give warnings when in-frame
    stop codons are found, invalid codons are found, or when the sequence length
    is not a multiple of three. Will translate the full sequence or until an
    unspecified or in-frame codon is encountered.

    :param sequence_dict: dictionary with gene IDs and CDS sequences
    :param skip_invalid: bool, skip invalid CDS? (default translates to first
        stop codon or end)
    :return: dictionary with gene IDs and proteins sequences
    """
    # TODO I should just use the Biopython translator
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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '',  'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '',  'TGG': 'W',
    }
    protein_dict = {}

    with ProgressBar(
            max_value=len(sequence_dict.keys()) + 1, redirect_stdout=True,
            redirect_stderr=True
    ) as pb:
        j = 0
        total = 0
        for key, val in sequence_dict.items():
            j += 1
            aa_seq = ''
            if len(val) % 3 != 0:
                print('Sequence length != multiple of 3 for {}!'.format(key))
                total += 1
            invalid = False
            for i in range(0, len(val), 3):
                if val[i:i + 3] not in aa_dict.keys():
                    print('Invalid codon {0:>3} in {1}'.format(val[i:i+3], key))
                    invalid = True
                    total += 1
                    break
                else:
                    if aa_dict[val[i:i + 3]] == '' and i+3 != len(val):
                        print('In-frame STOP codon in {0} at position {1}:{2}'
                              ''.format(key, i, i+3))
                        invalid = True
                        total += 1
                        break
                    aa_seq += aa_dict[val[i:i + 3]]
            if invalid and skip_invalid:
                continue
            protein_dict[key] = aa_seq
            pb.update(j)
    logging.warning("There were {} warnings during translation".format(total))
    return protein_dict


def write_fasta(seq_dict, output_file):
    """
    Write a sequence dictionary to a fasta file.

    :param seq_dict: sequence dictionary, see :py:func:`read_fasta`
    :param output_file: output file name
    """
    with open(output_file, 'w') as o:
        for key, val in seq_dict.items():
            o.write('>' + key + '\n')
            o.write(val + '\n')
    return output_file


def filter_one_vs_one_families(gene_families, s1, s2):
    """
    Filter one-vs-one ortholog containing families for two given species.

    :param gene_families: gene families fil in raw MCL format, see
        :py:func:`process_gene_families`
    :param s1: species 1 prefix
    :param s2: species 2 prefix
    :return: one-vs.-one ortholog containing gene families.
    """
    to_delete = []
    for key, val in gene_families.items():
        count = 0
        for gene in val:
            prefix = gene.split('|')[0]
            if prefix == s1 or prefix == s2:
                count += 1
        if count != 2:
            to_delete.append(key)
    for k in to_delete:
        del gene_families[k]
    return gene_families


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
        :param id_string: keyword for retrieving the gene ID from the 9th column
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
                    gene_dict = {x.split('=')[0]: x.split('=')[1] for x in
                                 gene_l if len(x.split('=')) == 2}

                    if chromosome not in self.genome:
                        self.genome[chromosome] = {}
                        self.gene_lists[chromosome] = []
                        self.colors[chromosome] = _random_color()

                    self.genome[chromosome][gene_dict[id_string]] = {
                        'orientation': orientation, 'start': start,
                        'stop': stop}
                    self.gene_lists[chromosome].append(
                            (gene_dict[id_string], orientation, start, stop))
        return

    def karyotype_json(self, out_file='genome.json'):
        """
        Generate karyotype data file in json format (as per Circos.js/d3.js)

        :param out_file: output file name
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


class gaussian_kde(object):
    """
    from: https://stackoverflow.com/questions/27623919/weighted-gaussian-kernel-
    density-estimation-in-python

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
        # http://surveyanalysis.org/wiki/Design_Effects_and_Effective_Sample_
        # Size#Kish.27s_approximate_formula_for_computing_effective_sample_size
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
                msg = "points have dimension %s, dataset has dimension %s" % \
                      (d, self.d)
                raise ValueError(msg)

        # compute the normalised residuals
        chi2 = cdist(points.T, self.dataset.T, 'mahalanobis',
                     VI=self.inv_cov) ** 2
        # compute the pdf
        result = np.sum(np.exp(-.5 * chi2) * self.weights,
                        axis=1) / self._norm_factor

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
        elif np.isscalar(bw_method) and not isinstance(bw_method, str):
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
            self._data_covariance = np.atleast_2d(
                    np.dot(_residual * self.weights, _residual.T))
            # Correct for bias (http://en.wikipedia.org/wiki/Weighted_arithmetic
            # _mean#Weighted_sample_covariance)
            self._data_covariance /= (1 - np.sum(self.weights ** 2))
            self._data_inv_cov = np.linalg.inv(self._data_covariance)

        self.covariance = self._data_covariance * self.factor ** 2
        self.inv_cov = self._data_inv_cov / self.factor ** 2
        self._norm_factor = np.sqrt(
                np.linalg.det(2 * np.pi * self.covariance))  # * self.n
