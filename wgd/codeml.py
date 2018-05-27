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

A python wrapper for codeml (PAML package, Yang 2007)

--------------------------------------------------------------------------------
"""
# TODO have a look at the BioPython wrapper

import pandas as pd
import numpy as np
import os
import subprocess
import re
import logging


def _write_control(f, control_dict):
    """
    Write the control file

    :param f: file handle
    :param control_dict: control keywords
    :return: written control file
    """
    for key in control_dict.keys():
        f.write('{0} = {1}\n'.format(key, control_dict[key]))


def _parse_codeml_out(codeml_out):
    """
    Parse the rst file to get pairwise Ks estimates

    :param codeml_out:
    :return: dictionary with CODEML results
    """
    if codeml_out is None:
        return None, None

    if not os.path.isfile(codeml_out):
        return None, None

    # compile regular expressions
    gene_pair_p = re.compile('\((.+)\) \.\.\. \d+ \((.+)\)')
    ks_p = re.compile('\s+dS\s*=\s*(\d+\.\d+)')
    ka_p = re.compile('\s+dN\s*=\s*(\d+\.\d+)')
    w_p = re.compile('\s+dN\/dS\s*=\s*(\d+\.\d+)')
    likelihood = re.compile('lnL\s*=(\D*\d+.\d+)')

    # read codeml output file
    with open(codeml_out, 'r') as f:
        file_content = f.read()
    codeml_results = file_content.split(
            'pairwise comparison')[-1].split("\n\n\n")[1:]

    columns = set()
    for pairwise_estimate in codeml_results:
        gene_1, gene_2 = gene_pair_p.search(pairwise_estimate).group(1), \
                         gene_pair_p.search(pairwise_estimate).group(2)
        columns.add(gene_1)
        columns.add(gene_2)

    results_dict = {
        'Ks': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns))),
        'Ka': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns))),
        'Omega': pd.DataFrame(
                np.zeros((len(list(columns)), len(list(columns)))),
                index=sorted(list(columns)),
                columns=sorted(list(columns)))
    }

    # populate results
    ln_l = None

    for pairwise_estimate in codeml_results:
        gene_1, gene_2 = gene_pair_p.search(pairwise_estimate).group(
                1), gene_pair_p.search(pairwise_estimate).group(2)
        ks_value_m = ks_p.search(pairwise_estimate)
        ka_value_m = ka_p.search(pairwise_estimate)
        w_m = w_p.search(pairwise_estimate)
        likelihood_m = likelihood.search(pairwise_estimate)
        if likelihood_m:
            if not ln_l:
                ln_l = 0
            try:
                ln_l += float(likelihood_m.group(1))
            except ValueError:
                logging.debug('Not a valid ln(L) value')
        else:
            logging.warning('No ln(L) value found!')

        # On the PLAZA 4.0 Vitis vinifera genome I had an issue with a pattern
        # match that was not retrieved. So now I check whether there is a match
        # and give a warning. anyway this shouldn't be necessary! if there is
        # codeml output, it should have the ks, ka and w values!
        if ks_value_m:
            ks_value = ks_value_m.group(1)
        else:
            logging.warning("No Ks value found in codeml file!")
            return None, None

        if ka_value_m:
            ka_value = ka_value_m.group(1)
        else:
            return None, None

        if w_m:
            w = w_m.group(1)
        else:
            return None, None

        results_dict['Ks'][gene_1][gene_2] = ks_value
        results_dict['Ks'][gene_2][gene_1] = ks_value
        results_dict['Ka'][gene_1][gene_2] = ka_value
        results_dict['Ka'][gene_2][gene_1] = ka_value
        results_dict['Omega'][gene_1][gene_2] = w
        results_dict['Omega'][gene_2][gene_1] = w

    return {'results': results_dict, 'raw': file_content}, ln_l


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

    Usage examples:

    Run codeml with default (Ks analysis) settings on a multiple sequence
    alignment (msa.fasta):::

        >>> codeml = Codeml()
        >>> codeml.run_codeml('msa.fasta')

    Change control setting CodonFreq to the F1x4 option (1) and rerun above
    example:::

        >>> codeml.control['CodonFreq'] = 1
        >>> codeml.run_codeml('msa.fasta')

    Do the above example directly without modifying the control settings in the
    dict directly:::

        >>> codeml = Codeml(CodonFreq=1)
        >>> codeml.run_codeml('msa.fasta')

    Print the default control settings:::

        >>> codeml = Codeml()
        >>> print(codeml)
    """

    def __init__(self, codeml='codeml', tmp='./', id='tmp', out_file=None,
                 **kwargs):
        """
        Codeml wrapper init. Initializes the default control file for Ks
        analysis as proposed by Vanneste et al. (2013). Takes as keyword
        arguments the options from the normal codeml distribution. Control
        settings are stored in a dictionary that can be accessed with the
        `.control` attribute

        :param codeml: path to codeml executable (by default will look for
            codeml in the system PATH)
        :param tmp: path to temporary directory, will default to the current
            working directory ('./')
        :param id: filename prefix for output/tmp files
        :param kwargs: any codeml control option (see PAML user guide)
        """
        self.codeml = codeml
        self.tmp = tmp

        if not os.path.isdir(self.tmp):
            raise NotADirectoryError(
                    'tmp directory {} not found!'.format(self.tmp))

        self.id = id
        self.control_file = os.path.join(self.tmp, self.id + '.ctrl')
        if not out_file:
            self.out_file = os.path.join(self.tmp, self.id + '.codeml')
        else:
            self.out_file = out_file
        self.control = {
            'seqfile': None, 'outfile': self.out_file,
            'noisy': 0, 'verbose': 0, 'runmode': -2, 'seqtype': 1,
            'CodonFreq': 2, 'clock': 0, 'aaDist': 0,
            'aaRatefile': 'dat/jones.dat', 'model': 0, 'NSsites': 0, 'icode': 0,
            'Mgene': 0, 'fix_kappa': 0,
            'kappa': 2, 'fix_omega': 0, 'omega': .4, 'fix_alpha': 1, 'alpha': 0,
            'Malpha': 0, 'ncatG': 8,
            'getSE': 0, 'RateAncestor': 1, 'Small_Diff': .5e-6, 'cleandata': 1,
            'method': 0
        }

        # update the control with kwargs
        for x in kwargs.keys():
            if x not in self.control:
                raise KeyError("{} is not a valid keyword for the codeml "
                               "control file.".format(x))
            else:
                self.control[x] = kwargs[x]

    def __str__(self):
        """
        String method for Codeml wrapper, prints current control settings

        :return: string representation of the control file
        """
        string = "codeml control 'file'\n---------------------\n"
        for key in sorted(self.control.keys()):
            string += '{0} = {1}\n'.format(key, self.control[key])
        return string

    def run_codeml(self, msa=None, raw=False, preserve=False, times=1):
        """
        Run codeml on a multiple sequence alignment file

        :param msa: multiple sequence alignment file
        :param raw: boolean, return raw results?
        :param preserve: boolean, preserve intermediate files?
        :param times: integer, perform codeml multiple times (average results)
        :return: dictionary with Ks, Kn and Kn/Ks (omega) values
        """
        if msa is not None:
            self.control['seqfile'] = msa

        elif msa is None and self.control['seqfile'] is None:
            raise ValueError("No sequence file provided!")

        # write the control file
        with open(self.control_file, 'w') as f:
            _write_control(f, self.control)

        # perform actual codeml times times
        output = []
        logging.debug("Performing codeml {} times".format(times))

        best = None
        best_index = 0
        for i in range(times):
            logging.debug(
                    "Codeml iteration {0} for {1}".format(str(i + 1), msa))

            # codeml hangs when filename is too long apparently
            # however in most use cases the files are in the working directory
            if os.path.isfile(os.path.basename(self.control_file)):
                self.control_file = os.path.basename(self.control_file)

            subprocess.run([self.codeml, self.control_file],
                           stdout=subprocess.PIPE)
            subprocess.run(
                    ['rm', '2ML.dN', '2ML.dS', '2ML.t', '2NG.dN', '2NG.dS',
                     '2NG.t', 'rst', 'rst1', 'rub'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if not os.path.isfile(self.out_file):
                logging.warning(
                        'Codeml output file {} not found'.format(self.out_file))
                return None

            d, likelihood = _parse_codeml_out(self.out_file)
            output.append(d)
            if not best or likelihood > best:
                best = likelihood
                best_index = i

        logging.debug('Best MLE: ln(L) = {}'.format(best))
        results = output[best_index]

        if not output:
            return None

        os.remove(self.control_file)

        if raw:
            return results

        else:
            if results is not None:
                return results['results'], os.path.abspath(self.out_file)
            return None, None



