import pandas as pd
import numpy as np
import subprocess as sp
import logging
import os
import re


def _write_aln_codeml(aln, fname):
    with open(fname, "w") as f:
        f.write("{} {}\n".format(len(aln), aln.get_alignment_length()))
        for s in aln:
            f.write("{}\n".format(s.id))
            f.write("{}\n".format(s.seq))

def _parse_pairwise(codeml_out):
    with open(codeml_out, "r") as f:
        content = f.read()
    x = content.split("pairwise comparison")[2].split("\n")
    y = [x for x in map(lambda i: _parse_pair(x, i), range(3, len(x), 7))]
    d = pd.DataFrame.from_dict(y).set_index("pair")
    return d

def _parse_pair(lines, start):
    x = lines[start].split()
    a = x[1].strip("()")
    b = x[4].strip("()")
    l = float(lines[start+1].split("=")[1].strip())
    x = lines[start+4].split("=")
    y = [x[0]]
    for elem in x[1:]:
        y += elem.split()
    x = {y[i].strip(): float(y[i+1].strip()) for i in range(0, len(y)-1, 2)}
    p = "__".join(sorted([a, b]))
    x.update({"pair": p, "gene1": a, "gene2": b, "l": l})
    return x


class Codeml:
    """
    Class for codeml (PAML Yang 2007) python wrapper. Defines the controle file
    and enables running codeml from within python in one line of code.

    :param codeml: path to codeml executable (by default will look for codeml
        in the system PATH)
    :param tmp: path to temporary directory, will default to the current working
        directory ('./')
    :param id: filename prefix for output/tmp files
    :param kwargs: other codeml keyword arguments (see PAML user guide)::

        'seqfile': None,
        'outfile': self.out_file,
        'noisy': 0,
        'verbose': 0,
        'runmode': -2,
        'seqtype': 1,
        'CodonFreq': 2,
        'clock': 0,
        'aaDist': 0,
        'aaRatefile': 'dat/jones.dat',
        'model': 0,
        'NSsites': 0,
        'icode': 0,
        'Mgene': 0,
        'fix_kappa': 0,
        'kappa': 2,
        'fix_omega': 0,
        'omega': .4,
        'fix_alpha': 1,
        'alpha': 0,
        'Malpha': 0,
        'ncatG': 8,
        'getSE': 0,
        'RateAncestor': 1,
        'Small_Diff': .5e-6,
        'cleandata': 1,
        'method': 0
    """
    def __init__(self, aln, exe='codeml', tmp='./', prefix='codeml', **kwargs):
        """
        Codeml wrapper init. Initializes the default control file for Ks
        analysis as proposed by Vanneste et al. (2013). Takes as keyword
        arguments the options from the normal codeml distribution. Control
        settings are stored in a dictionary that can be accessed with the
        `.control` attribute

        :param aln: alignment file
        :param exe: path to codeml executable (by default will look for
            codeml in the system PATH)
        :param tmp: path to temporary directory, will default to the current
            working directory ('./')
        :param prefix: filename prefix for output/tmp files
        :param kwargs: any codeml control option (see PAML user guide)
        """
        if not os.path.isdir(tmp):
            raise NotADirectoryError('tmp directory {} not found!'.format(tmp))
        self.prefix = prefix
        self.aln = aln
        self.exe = exe
        self.tmp = tmp
        self.control_file = self.prefix + '.ctrl'
        self.aln_file = self.prefix + '.cdsaln'
        self.out_file = self.prefix + '.codeml'
        self.results = None
        self.control = {
            'seqfile': self.aln_file, 
            'outfile': self.out_file,
            'noisy': 0, 
            'verbose': 0, 
            'runmode': -2, 
            'seqtype': 1,
            'CodonFreq': 2, 
            'clock': 0, 
            'aaDist': 0,
            'aaRatefile': 'dat/jones.dat', 
            'model': 0, 
            'NSsites': 0, 
            'icode': 0,
            'Mgene': 0, 
            'fix_kappa': 0,
            'kappa': 2, 
            'fix_omega': 0, 
            'omega': .4, 
            'fix_alpha': 1, 
            'alpha': 0,
            'Malpha': 0, 
            'ncatG': 8,
            'getSE': 0, 
            'RateAncestor': 1, 
            'Small_Diff': .5e-6, 
            'cleandata': 1,
            'method': 0}
        # update the control with kwargs
        for x in kwargs.keys():
            if x not in self.control:
                raise KeyError("{} is not a valid codeml param.".format(x))
            else:
                self.control[x] = kwargs[x]

    def __str__(self):
        """
        String method for Codeml wrapper, prints current control settings

        :return: string representation of the control file
        """
        x = ['{0} = {1}\n'.format(k, v) for (k,v) in sorted(self.control.items())]
        return "\n".join(x)

    def write_ctrl(self):
        with open(self.control_file, "w") as f:
            f.write(str(self))

    def run_codeml(self, preserve=False, times=1):
        """
        Run codeml on a multiple sequence alignment file

        :param preserve: boolean, preserve intermediate files?
        :param times: integer, perform codeml multiple times (average results)
        :return: dictionary with Ks, Kn and Kn/Ks (omega) values
        """
        parentdir = os.path.abspath(os.curdir)
        os.chdir(self.tmp)
        with open(self.control_file, 'w') as f:
            self.write_ctrl()
        _write_aln_codeml(self.aln, self.aln_file)
        logging.debug("Performing codeml {} times".format(times))
        
        max_results = None 
        max_likelihood = None
        for i in range(times):
            logging.debug("Codeml iteration {0} for {1}".format(str(i+1), self.prefix))
            sp.run([self.exe, self.control_file], stdout=sp.PIPE)
            sp.run(['rm', '2ML.dN', '2ML.dS', '2ML.t', '2NG.dN', '2NG.dS',
                '2NG.t', 'rst', 'rst1', 'rub'], stdout=sp.PIPE, stderr=sp.PIPE)
            if not os.path.isfile(self.out_file):
                raise FileNotFoundError('Codeml output file not found')
            results = _parse_pairwise(self.out_file)
            likelihood = sum(results["l"])
            if not max_likelihood or likelihood > max_likelihood:
                max_likelihood = likelihood
                max_results = results

        logging.debug('Best MLE: ln(L) = {}'.format(max_likelihood))
        os.remove(self.control_file)
        if not preserve:
            os.remove(self.out_file)
        os.chdir(parentdir)
        self.results = max_results
        return max_results
