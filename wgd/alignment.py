#!/usr/bin/python3.5
"""
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

Alignment related tools.

Contribution note: It is trivial to add support for another aligner (besides
MUSCLE and PRANK). Especially when it's for amino acid alignments (codon
alignment as PRANK are a bit more involved). To add another aligner, only the
:py:meth:`MSA.run_aligner` method should be modified to include the aligner of
interest.
"""
# TODO: use BioPython alignIO stuff + rewrite stats for alignments

from .utils import read_fasta, write_fasta
import os
import subprocess
import logging
import itertools


# REWRITE
def prepare_aln(msa_file, nuc_seqs):
    """
    Main wrapper function for alignments for Ks distributions, takes an 
    alignment file as input (protein or nucleotide), and returns a nucleotide 
    alignment (codon) file and pairwise alignment statistics.

    :param msa_file: multiple sequence alignment file path
    :param nuc_seqs: nucleotide sequences (set to ``None`` if input is a 
        nucleotide alignment).
    :return: file path to nucleotide alignment, pairwise alignment statistics
    """
    with open(msa_file, 'r') as f:
        aln = read_fasta(msa_file)
        if not nuc_seqs:
            aln = pal2nal(aln, nuc_seqs)
        stats = pairwise_alignment_stats_(aln)
    out_path = msa_file + '.nuc'
    write_alignment_codeml(aln, out_path)
    return out_path, stats


def pal2nal(pal, nuc_seqs):
    """
    Protein alignment to nucleotide alignment converter.

    :param pal: protein alignment dictionary
    :param nuc_seqs: nucleotide sequences
    :return: nucleotide alignment dictionary
    """
    nal = {}
    nal_ = ''
    for pid, seq in pal.items():
        if pid not in nuc_seqs:
            logging.warning(
                "Sequence {} in protein alignment not found in CDS sequences"
                "".format(pid)
            )
        else:
            nuc = nuc_seqs[pid]
            nal_ = ''
            j = 0
            for i in range(len(seq)):
                if seq[i] == '-':
                    nal_ += '---'
                else:
                    nal_ += nuc[j:j + 3]
                    j += 3
            nal[pid] = nal_
    return nal


def pairwise_alignment_stats_(aln):
    """
    Get pairwise alignment statistics.

    :param aln: alignment dictionary
    :return: dictionary with pairwise statistics
    """
    pairs = itertools.combinations(list(aln.keys()), 2)
    stats_dict = {}
    aln_len = len(list(aln.values())[0])

    # loop over all pairs in the alignment
    for pair in pairs:
        id1, id2 = pair
        s1, s2 = strip_gaps_pair(aln[id1], aln[id2])
        identity = (len(s1) - hamming_distance(s1, s2)) / len(s1)
        stats_dict["_".join(sorted(pair))] = {
            "AlignmentIdentity": identity,
            "AlignmentLength": aln_len,
            "AlignmentLengthStripped": len(s1),
            "AlignmentCoverage": len(s1)/aln_len
        }
    return stats_dict


def strip_gaps_pair(s1, s2):
    """
    Strip gaps for an aligned sequence pair.

    :param s1: sequence 1
    :param s2: sequence 2
    :return: two stripped sequences
    """
    s1_, s2_ = '', ''
    for i in range(len(s1)):
        if s1[i] == '-' or s2[i] == '-':
            continue
        else:
            s1_ += s1[i]
            s2_ += s2[i]
    return s1_, s2_


def hamming_distance(s1, s2):
    """
    Return the Hamming distance between equal-length sequences
    
    :param s1: string 1
    :param s2: string 2
    :return: the Hamming distances between s1 and s2
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


class MSA:
    """
    Multiple multiple sequence alignment (MSA) programs python wrappers.
    Currently only runs with default settings.
    """

    def __init__(self, muscle='muscle', prank='prank', tmp='./'):
        """
        Muscle wrapper init.

        :param muscle: path to Muscle executable, will by defult look for Muscle
            in the system PATH
        :param prank: path to PRANK executable
        :param tmp: directory to store temporary files.
        """
        self.muscle = muscle
        self.prank = prank
        self.tmp = tmp

        if not os.path.isdir(self.tmp):
            raise NotADirectoryError(
                    'tmp directory {} not found!'.format(self.tmp))

    def run_aligner(self, sequences, aligner='muscle', file=None,
                    return_dict=False):
        """
        Run MSA on sequences stored in dictionary with a particular aligner

        :param sequences: dictionary of sequences or input fasta file
        :param file: output file name, if not provided output is dictionary
        :param aligner: alignment software to use (prank|muscle)
        :return: MSA file
        """
        if not file:
            file_name = 'msa'
        else:
            file_name = file

        # sequences provided in dictionary
        if type(sequences) == dict:
            target_path_fasta = os.path.join(
                self.tmp, '{}.fasta'.format(file_name)
            )
            write_fasta(sequences, target_path_fasta)

        # sequences provided as fasta file
        elif os.path.isfile(sequences):
            target_path_fasta = sequences

        # error
        else:
            logging.error('{} is dictionary nor file path'.format(sequences))
            return None

        target_path_msa = os.path.join(self.tmp, '{}.msa'.format(file_name))

        if aligner == 'muscle':
            subprocess.run([
                self.muscle, '-quiet', '-in', target_path_fasta, '-out',
                target_path_msa
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        elif aligner == 'prank':
            subprocess.run([self.prank, '-codon', '-d=' + target_path_fasta,
                            '-o=' + target_path_msa],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subprocess.run(['mv', '{}.best.fas'.format(target_path_msa),
                            target_path_msa])

        subprocess.run(['rm', target_path_fasta], stdout=subprocess.PIPE)

        if return_dict:
            out = read_fasta(target_path_msa)
            return target_path_msa, out

        return target_path_msa
