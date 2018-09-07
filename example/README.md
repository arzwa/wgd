# Examples

This directory is designated for example analyses and data. You can find
the _Arabidopsis thaliana_ Ks distribution for the full paranome as well
as anchors in the data directory. These are used in the example
notebook(s).

The example data files were obtained using the commands available in the
wgd CLI. The input data was retrieved from the PLAZA (4.0) platform. To
obtain these results, see below.

## Getting the results yourself

Get the sequence data

    wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04//Fasta/cds.ath.fasta.gz
    gunzip cds.ath.fasta.gz

Note that this annotation contains mitochondrial and chloroplast genes
as well as transposable elements. To skip these do

     grep -A 1 ">AT[1-9]G" cds.ath.fasta > cds.ath_filtered.fasta

For a quick test analysis, get a random sample of sequences from the
file (skip this if you want to do the full analysis).

    grep ">" cds.ath_filtered.fasta | shuf | head -n 1000 > ids
    grep -A 1 -f ids cds.ath_filtered.fasta | sed "/--/d" > sample.fasta
    rm ids

Run an all-_vs._-all Blastp analysis and cluster using MCL, assuming
you're working with `sample.fasta` (replace with `cds.ath_filtered.fasta`
for a full analysis).

    wgd mcl -s sample.fasta --cds --mcl

A directory named `wgd_blast` will appear that will contain some gene
families in the file `sample.fasta.blast.tsv.mcl`. Let's move the file
to the work dir and give it a shorter name

    mv wgd_blast/sample.fasta.blast.tsv.mcl ./sample.mcl

Now let's compute a Ks distribution (use `-n` to set the number of cores
to use, defaults to 4).

    wgd ksd sample.mcl sample.fasta

This creates a directory `wgd_ksd`. You can check the histograms that
were generated:

    display wgd_ksd/*svg

or inspect the Ks distribution itself:

    head -n 3 wgd_ksd/*tsv
```
    AlignmentCoverage	AlignmentIdentity	AlignmentLength	AlignmentLengthStripped	Distance	Family	Ka	Ks	Node	Omega	Outlier	Paralog1	Paralog2	WeightOutliersExcluded	WeightOutliersIncluded
AT1G77815__AT3G09510	0.18941	0.50538	1473.0	279.0	1.31918	GF_000060	0.7452	2.4193	2.0	0.308	False	AT1G77815	AT3G09510	1.0	1.0
AT2G27980__AT2G37520	0.73607	0.49256	3285.0	2418.0	1.21274	GF_000074	0.6217	10.1451	2.0	0.0613	True	AT2G27980	AT2G37520	0.0	1.0
```

To get anchor pairs, first download the GFF file

    wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04//GFF/ath/Arabidopsis_thaliana.COL0.Araport11.longest_transcript.all_features.gff3.gz
    gunzip Arabidopsis_thaliana.COL0.Araport11.longest_transcript.all_features.gff3.gz
    mv Arabidopsis_thaliana.COL0.Araport11.longest_transcript.all_features.gff3 ath.gff

Then run `wgd syn`:

    wgd syn -ks wgd_ksd/sample.fasta.ks.tsv -f gene -a ID ath.gff sample.mcl

Chances are huge that you get the `WARNING	No multiplicons found!`
warning when using the small test data, in that case nothing interesting
happens. However when doing the full analysis, you should find dotplots
and anchor-pair Ks distributions in the `wgd_syn` directory.


