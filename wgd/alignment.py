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
# TODO: use BioPython alignIO stuff

from .utils import read_fasta
import os
import subprocess
import logging
import itertools


def _nucleotide_msa_from_protein_msa(
        protein_msa_sequences, nucleotide_sequences):
    """
    Make a nucleotide multiple sequence alignment from a protein multiple
    sequence alignment

    :param protein_msa_sequences: dictionary with protein sequences (aligned)
    :param nucleotide_sequences: dictionary with nucleotide sequences
        (unaligned)
    :return: dictionary with nucleotide sequences (aligned)
    """
    nucleotide_msa = {}
    lengths = set()

    for seq in protein_msa_sequences.values():
        length = len(seq)
        lengths.add(length)

    if len(list(lengths)) != 1:
        logging.warning(
            "Not all sequences have the same length after alignment!")
        return None

    for gene_ID in protein_msa_sequences.keys():
        protein = protein_msa_sequences[gene_ID]
        if gene_ID not in nucleotide_sequences.keys():
            logging.warning(
                "Gene {} not found in nucleotide fasta!".format(gene_ID))
        else:
            nucleotide = nucleotide_sequences[gene_ID]
            nucleotide_msa_sequence = ''
            j = 0
            for i in range(len(protein)):
                if protein[i] == '-':
                    nucleotide_msa_sequence += '---'
                else:
                    nucleotide_msa_sequence += nucleotide[j:j + 3]
                    j += 3
            nucleotide_msa[gene_ID] = nucleotide_msa_sequence

    return nucleotide_msa


def _strip_gaps(msa_dict):
    """
    Strip gap positions from a multiple sequence alignment (MSA).
    Finds the positions in the strings that have a gap and removes them in all
    sequences.

    :param msa_dict: dictionary with aligned nucleotide sequences.
    """
    if msa_dict is None:
        return None
    if len(list(msa_dict.values())) == 0:
        return None

    indices = set()
    length = len(list(msa_dict.values())[0])
    for i in range(length):
        for gene_ID in msa_dict.keys():
            # Prevent index error, raise warning instead
            if i >= len(msa_dict[gene_ID]):
                logging.warning(
                    "Seems there are unequal string lengths after alignment in "
                    "this gene family. Occurred at gene: {}.".format(gene_ID))
                return None
            if msa_dict[gene_ID][i] == '-':
                indices.add(i)
    indices = list(indices)

    for gene_ID in msa_dict.keys():
        sequence_list = [i for i in msa_dict[gene_ID]]
        for index in sorted(indices, reverse=True):
            del sequence_list[index]
        msa_dict[gene_ID] = ''.join(sequence_list)

    return msa_dict


def multiple_sequence_aligment_nucleotide(msa_protein, nucleotide_sequences,
                                          min_length=100, aligner='muscle'):
    """
    Make a nucleotide multiple sequence alignment based on a protein alignment

    :param msa_protein: dictionary of aligned protein sequences
    :param nucleotide_sequences: nucleotide sequence dictionary
    :param min_length: minimum alignment length to consider
    :param aligner: alignment program used, if prank, than msa_protein will be
        interpreted as codon alignment
    :return: nucleotide MSA, length, stripped length
    """
    if not os.path.isfile(msa_protein):
        logging.warning('MSA file {} not found!'.format(msa_protein))
        return None

    # Back-translate -----------------------------------------------------------
    # if aligner is prank than the family is codon aligned so no need to
    # back-translate
    protein_msa_sequences = read_fasta(msa_protein)
    if aligner != 'prank':
        nucleotide_msa = _nucleotide_msa_from_protein_msa(protein_msa_sequences,
                                                          nucleotide_sequences)
    else:
        nucleotide_msa = protein_msa_sequences

    if nucleotide_msa is None:
        return None

    # get alignment statistics for quality control -----------------------------
    alignment_stats = pairwise_alignment_stats(nucleotide_msa)
    l_unstripped = len(list(nucleotide_msa.values())[0])

    # strip gaps ---------------------------------------------------------------
    nucleotide_msa = _strip_gaps(nucleotide_msa)

    if nucleotide_msa is None:  # probably not necessary
        return None

    elif len(list(nucleotide_msa.values())[0]) < min_length:
        return None

    l_stripped = len(list(nucleotide_msa.values())[0])

    msa_nuc = msa_protein + '.nuc'
    with open(msa_nuc, 'w') as o:
        o.write("\t{0}\t{1}\n".format(len(nucleotide_msa.keys()),
                                      len(list(nucleotide_msa.values())[0])))
        for gene_ID in nucleotide_msa.keys():
            o.write("{}\n".format(gene_ID))
            o.write(nucleotide_msa[gene_ID])
            o.write("\n")

    return [msa_nuc, int(l_unstripped), int(l_stripped), alignment_stats]


def get_pairwise_nucleotide_alignments(msa_protein, nucleotide_sequences,
                                       min_length=100, aligner='muscle'):
    """

    :param msa_protein:
    :param nucleotide_sequences:
    :return: a list with file names
    """
    if not os.path.isfile(msa_protein):
        logging.warning('MSA file {} not found!'.format(msa_protein))
        return None

    pairwise_alignments = []

    # back-translate -----------------------------------------------------------
    # if aligner is prank than the family is codon aligned so no need to
    # back-translate
    protein_msa_sequences = read_fasta(msa_protein)
    if aligner != 'prank':
        nucleotide_msa = _nucleotide_msa_from_protein_msa(protein_msa_sequences,
                                                          nucleotide_sequences)
    else:
        nucleotide_msa = protein_msa_sequences

    if nucleotide_msa is None:
        return None

    # get all alignment pairs with gap-stripped length > min_length ------------
    pairs = itertools.combinations(list(nucleotide_msa.keys()), 2)
    i = 1
    for pair in pairs:
        pair_dict = {pair[0]: nucleotide_msa[pair[0]],
                     pair[1]: nucleotide_msa[pair[1]]}
        # strip gaps
        l_unstripped = len(list(pair_dict.values())[0])
        pair_dict = _strip_gaps(pair_dict)

        if pair_dict is None:  # probably not necessary
            continue

        elif len(list(pair_dict.values())[0]) < min_length:
            logging.warning('Stripped alignment for gene {0} and gene {1} '
                            'too short (< {2} nucleotides)'.format(
                    pair[0], pair[1], min_length
            ))
            continue

        l_stripped = len(list(pair_dict.values())[0])
        alignment_stats = pairwise_alignment_stats(pair_dict, pair=True)

        # write alignment
        file_name = msa_protein + '.' + str(i) + '.nuc'
        write_alignment_codeml(pair_dict, file_name)
        pairwise_alignments.append((file_name, pair[0], pair[1],
                                    l_stripped, l_unstripped, alignment_stats))
        i += 1

    return pairwise_alignments


def write_alignment_codeml(alignment, file_name):
    with open(file_name, 'w') as o:
        o.write("\t{0}\t{1}\n".format(
                len(alignment.keys()), len(list(alignment.values())[0])))
        for gene_ID in alignment.keys():
            o.write("{}\n".format(gene_ID))
            o.write(alignment[gene_ID])
            o.write("\n")


def pairwise_alignment_stats(msa, pair=False):
    """
    Get pairwise stats from MSA.
    Return as  dictionary {Gene1: {Gene2: [0.8, 0.7], Gene3: [...], ...}

    :param msa: MSA dictionary
    :return: dictionary
    """
    stats = {}
    stats_pair = None
    seqs = [(x, msa[x]) for x in sorted(msa.keys())]

    for i in range(len(seqs)):
        stats[seqs[i][0]] = {}
        for j in range(i, len(seqs)):
            id1, id2 = seqs[i][0], seqs[j][0]
            d = _strip_gaps({id1: seqs[i][1], id2: seqs[j][1]})
            try:
                # % identity
                stats_pair = [(len(d[id1]) - hamming_distance(d[id1], d[id2]))
                              / len(d[id1]), len(d[id1]) / len(seqs[i][1])]
                # coverage
                stats[id1][id2] = stats_pair
            except:
                stats[id1][id2] = [0, 0]
    if pair:
        return stats_pair

    return stats


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


class MSA:
    """
    Multiple multiple sequence alignment (MSA) programs python wrappers.
    Currently only runs with default settings.

    Usage examples

    * Sequences in dictionary, output as dictionary::

        >>> MSA().run_aligner({'bear_gene': 'BEARBEAR', 'hare_gene': 'HAREERAH', 'yeast_gene': 'BEERBEARBEER'})
        {'bear_gene': '----BEARBEAR-',
        'hare_gene': '-----HAREERAH',
        'yeast_gene': 'BEERBEARBEER-'}

    * Sequences in fasta file, output as fasta file ``(msa.fasta)``::

        >>> msa = MSA()
        >>> msa.run_aligner('./sequences.fasta', './msa.fasta')

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
            target_path_fasta = os.path.join(self.tmp,
                                             '{}.fasta'.format(file_name))

            with open(target_path_fasta, 'w') as o:
                for gene in sequences.keys():
                    o.write('>{}\n'.format(gene))
                    o.write(sequences[gene])
                    o.write('\n')

        # sequences provided as fasta file
        elif os.path.isfile(sequences):
            target_path_fasta = sequences

        # error
        else:
            logging.error(
                '{} is not a dictionary and also not a file path'.format(
                    sequences))
            return None

        target_path_msa = os.path.join(self.tmp, '{}.msa'.format(file_name))

        if aligner == 'muscle':
            subprocess.run(
                    [self.muscle, '-quiet', '-in', target_path_fasta, '-out',
                     target_path_msa],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
