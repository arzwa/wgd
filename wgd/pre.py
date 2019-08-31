from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import logging

def check_cds(cdsfile, out1, out2, rename=False, prefix=None):
    """
    Check a file of CDS sequences whether they are proper CDS.

    :param cdsfile: CDS file name
    :param out1: output filename for proper CDS sequences
    :param out2: output filename for dubitable CDS sequences
    :param rename: should the gene IDs be renamed?
    :param prefix: prefix in case rename=True (defaults to the filename)
    :return: nada
    """
    if not prefix:
        prefix = cdsfile
    x = 0
    y = 0
    d = {}
    with open(out2, "w") as f2:
        with open(out1, "w") as f1:
            for i, seq_record in enumerate(SeqIO.parse(cdsfile, 'fasta')):
                if rename:
                    gid = "{0}_{1:0>5}".format(prefix, i)
                    d[seq_record.id] = gid
                else:
                    gid = seq_record.id
                try:
                    aa_seq = seq_record.seq.translate(to_stop=True, cds=True)
                except TranslationError as e:
                    logging.error("Translation error ({}) in sequence "
                                  "{}".format(e, seq_record.id))
                    f2.write(">{}\n{}\n".format(gid, seq_record.seq))
                    y += 1
                    continue
                f1.write(">{}\n{}\n".format(gid, seq_record.seq))
                x += 1
    t = x + y
    logging.info("{}/{} ({:.2f}%) sequences are perfect CDS (in {})".format(
        x, t, 100*x/t, out1))
    logging.info("{}/{} ({:.2f}%) sequences are not perfect CDS (in {})".format(
        y, t, 100*y/t, out2))
    if rename:
        with open(cdsfile + ".ids.csv", "w") as f:
            f.write("{},{}\n".format(cdsfile, out1))
            for k,v in d.items():
                f.write("{},{}\n".format(k, v))
