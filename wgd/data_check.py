"""
`wgd` WGD analysis in python

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
This is a module for performing common data checks for paranoids.
--------------------------------------------------------------------------------
"""
from Bio import SeqIO, Seq
import logging


def check_cds(cds_file):
    """"""
    seq_dict = SeqIO.index(cds_file, 'fasta')
    invalid = []
    for k, v in seq_dict.items():
        try:
            v.seq.translate(cds=True)
        except Seq.CodonTable.TranslationError as e:
            logging.warning('{} is not a valid coding sequence: '.format(k))
            logging.warning(e)
            invalid.append(v)
    logging.info('Invalid CDS: {0:>6}/{1:>6} ({2:>4.2f}%)'.format(
            len(invalid), len(seq_dict), len(invalid)/len(seq_dict) * 100))
    SeqIO.write(invalid, cds_file + '.invalid_cds', 'fasta')