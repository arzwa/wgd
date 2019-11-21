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
"""
from .utils import read_fasta, write_fasta, log_subprocess
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
        if nuc_seqs:
            aln = pal2nal(aln, nuc_seqs)
        stats = pairwise_alignment_stats(aln)
    out_path = msa_file + '.nuc'
    success = write_alignment_codeml(aln, out_path)
    return out_path, stats, success


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


def get_pairwise_alns(aln, nuc_seqs, min_length=3):
    """
    Get all pairwise alignments and pairwise statistics.

    :param aln: alignment file
    :param nuc_seqs: nucleotide sequences dictionary
    :param min_length: minimum stripped alignment length necessary to include
        a pair
    """
    aln = read_fasta(aln)
    pairwise_alns = []
    pairs = itertools.combinations(list(aln.keys()), 2)
    stats_dict = {}
    for pair in pairs:
        id1, id2 = pair
        pid = '__'.join(sorted([id1, id2]))
        seqs = {x: aln[x] for x in (id1, id2)}
        if nuc_seqs:
            seqs = pal2nal(seqs, nuc_seqs)
        stats_dict[pid] = get_stats(seqs[id1], seqs[id2])
        s1, s2 = strip_gaps_pair(seqs[id1], seqs[id2])
        if len(s1) < min_length or len(s2) < min_length:
            continue
        pairwise_alns.append((pid, {id1: s1, id2: s2}))
    return pairwise_alns, stats_dict


def pairwise_alignment_stats(aln):
    """
    Get pairwise alignment statistics.

    :param aln: alignment dictionary
    :return: dictionary with pairwise statistics
    """
    pairs = itertools.combinations(list(aln.keys()), 2)
    stats_dict = {}

    # loop over all pairs in the alignment
    for pair in pairs:
        id1, id2 = pair
        stats_dict['__'.join(sorted(pair))] = get_stats(aln[id1], aln[id2])
    return stats_dict


def get_stats(s1, s2):
    s1_, s2_ = strip_gaps_pair(s1, s2)
    if len(s1_) == 0:
        identity = 0
    else:
        identity = (len(s1_) - hamming_distance(s1_, s2_)) / len(s1_)
    return {
        "AlignmentIdentity": identity,
        "AlignmentLength": len(s1),
        "AlignmentLengthStripped": len(s1_),
        "AlignmentCoverage": len(s1_)/len(s1)
    }

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


def write_alignment_codeml(alignment, file_name):
    """
    Write an alignment file for codeml

    :param alignment: alignment dictionary
    :param file_name: output file name
    """
    n = len(alignment.keys())
    if n == 0:
        return False
    m = len(list(alignment.values())[0])
    with open(file_name, 'w') as o:
        o.write("\t{0}\t{1}\n".format(n, m))
        for gene_ID in alignment.keys():
            o.write("{}\n".format(gene_ID))
            o.write(alignment[gene_ID])
            o.write("\n")
    return True


def align_muscle(in_file, out_file, bin_path='muscle'):
    """
    Perform multiple sequence alignment with MUSCLE

    :param bin_path: path to MUSCLE executable
    :param in_file: input fasta file
    :param out_file: output fasta file
    """
    cmd = [bin_path, '-in', in_file, '-out', out_file]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log_subprocess("MUSCLE", out)
    return out_file


def align_prank(in_file, out_file, bin_path='prank'):
    """
    Perform multiple sequence alignment with PRANK

    :param bin_path: path to PRANK executable
    :param in_file: input fasta file
    :param out_file: output fasta file
    """
    cmd = [bin_path, '-codon', '-d=' + in_file, '-o=' + out_file]    
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log_subprocess("PRANK", out)
    return out_file


def align_mafft(in_file, out_file, bin_path='mafft'):
    """
    Perform multiple sequence alignment with MAFFT

    :param bin_path: path to MAFFT executable
    :param in_file: input fasta file
    :param out_file: output fasta file
    """
    cmd = [bin_path, "--amino", '--maxiterate', '1000', '--localpair', in_file]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log_subprocess("MAFFT", out)
    with open(out_file, 'w') as f: f.write(out.stdout.decode('utf-8'))
    return out_file


def align(in_file, out_file, aligner):
    """
    Alignment wrapper

    :param aligner: alignment program
    :param in_file: input fasta file
    :param out_file: output fasta file
    """
    if not os.path.isfile(in_file):
        logging.warning("File not found {}".format(in_file))
        return None
    if aligner == 'prank':
        out = align_prank(in_file, out_file)
    elif aligner == 'mafft':
        out = align_mafft(in_file, out_file)
    else:
        out = align_muscle(in_file, out_file)
    return out
