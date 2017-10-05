#!/usr/bin/python3.5
"""
Arthur Zwaenepoel

GFF parsing tools. Defines a class 'Genome' which is a generic 
container for structural annotations
"""
import random
import json


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

    def parse_plaza_gff(self, gff_file, keyword='mRNA'):
        """
        Parse a PLAZA annotation file into a genome dictionary

        :param gff_file: input gff (PLAZA style)
        :param keyword: keyword for elements to parse out
        """
        self.parent_file = gff_file

        with open(gff_file, 'r') as f:
            for line in f:
                line = line.strip().split('\t')

                if line[2] == keyword:
                    chromosome = line[0]
                    start = line[3]
                    stop = line[4]
                    orientation = line[6]
                    gene = line[8].split(';')[0].split('=')[1]

                    if chromosome not in self.genome:
                        self.genome[chromosome] = {}
                        self.gene_lists[chromosome] = []
                        self.colors[chromosome] = _random_color()

                    self.genome[chromosome][gene] = {
                        'orientation': orientation, 'start': start, 'stop': stop}
                    self.gene_lists[chromosome].append(
                        (gene, orientation, start, stop))
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
