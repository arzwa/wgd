#!/usr/bin/python3.5
"""
Arthur Zwaenepoel
"""
# TODO: This should be adapted such that multiple genomes can be compared

import os
import subprocess
import logging
import pandas as pd
import matplotlib.pyplot as plt


# WRITE FILES AND CONFIG

def write_gene_lists(genome, output_dir='gene_lists'):
    """
    Write out the gene lists

    :param genome: Genome object 
    :param output_dir: output directory
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # meanwhile, store all genes in a set to later add as singletons
    all_genes = set()

    for chromosome in genome.genome.keys():
        with open(os.path.join(output_dir, chromosome + '.lst'), 'w') as o:
            for gene in genome.gene_lists[chromosome]:
                o.write(gene[0] + gene[1] + '\n')
                all_genes.add(gene[0])

    return all_genes


def write_families_file(families, all_genes, output_file='families.tsv'):
    """
    Write out families file

    :param families: MCL output 
    """
    counter = 1
    genes_seen = set()

    with open(families, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                line = line.strip().split('\t')
                for gene in line:
                    o.write(gene + '\t' + str(counter) + '\n')
                    genes_seen.add(gene)
                counter += 1

            # add genes not seen to the families file as singletons
            # I-ADHoRe throws an error when there are genes in the gene
            # lists that are not in the blast table/families file
            rest = all_genes - genes_seen
            for gene in list(rest):
                o.write(gene + '\t' + str(counter) + '\n')
                counter += 1
    return


def write_config_adhore(gene_lists, families, config_file_name='i-adhore.conf',
                        genome='genome', output_path='i-adhore_out', gap_size=30,
                        cluster_gap=35, q_value=0.75, prob_cutoff=0.01, anchor_points=3,
                        alignment_method='gg2', level_2_only='false', table_type='family',
                        multiple_hypothesis_correction='FDR', visualizeGHM='false',
                        visualizeAlignment='true', **kwargs):
    """
    Write out the config file for I-ADHoRe. See I-ADHoRe manual for information on
    parameter settings.

    :param gene_lists: directory with gene lists per chromosome
    :param families: file with gene to family mapping
    :return: configuration file
    """
    with open(config_file_name, 'w') as o:
        o.write('genome= {}\n'.format(genome))

        counter = 1
        for l in sorted(os.listdir(gene_lists)):
            o.write("{0} {1}/{2}\n".format(l[:-4], gene_lists, l))
            counter += 1

        o.write('blast_table= {}\n'.format(families))
        o.write('output_path= {}\n'.format(output_path))
        o.write('gap_size= {}\n'.format(gap_size))
        o.write('q_value= {}\n'.format(q_value))
        o.write('cluster_gap= {}\n'.format(cluster_gap))
        o.write('prob_cutoff= {}\n'.format(prob_cutoff))
        o.write('anchor_points= {}\n'.format(anchor_points))
        o.write('alignment_method= {}\n'.format(alignment_method))
        o.write('level_2_only= {}\n'.format(level_2_only))
        o.write('table_type= {}\n'.format(table_type))
        o.write('multiple_hypothesis_correction= {}\n'.format(
            multiple_hypothesis_correction))
        o.write('visualizeGHM= {}\n'.format(visualizeGHM))
        #o.write('visGPairs= {}\n'.format(output_path))
        o.write('visualizeAlignment= {}\n'.format(visualizeAlignment))

    return


def get_anchor_pairs(anchors_file, ks_file=None, out_file='anchors_ks.csv', species='SOI'):
    """
    Get anchor pairs and their corresponding Ks values (if provided)

    :param anchors_file: anchorpoints.txt output from I-ADHoRe 3.0
    :param ks_file: Ks calculations file from wgd toolkit (all.csv)
    :return: pandas dataframe(s): anchors and data frame
    """
    if not os.path.exists(anchors_file):
        logging.error('Anchor points file: `{}` not found'.format(anchors_file))
    else:
        logging.info('Anchor points file found.')

    anchors = pd.read_csv(anchors_file, sep='\t', index_col=0)[['gene_x', 'gene_y']]

    # give the pairs an identifier
    ids = []
    for x in anchors.index:
        ids.append("-".join(sorted([anchors.loc[x]['gene_x'], anchors.loc[x]['gene_y']])))
    anchors['pair_id'] = ids

    if ks_file:
        if not os.path.exists(ks_file):
            logging.error('Ks values file: `{}` not found'.format(ks_file))
        else:
            logging.info('Ks values file found.')

        ks = pd.read_csv(ks_file, index_col=0)

        # reindex by pair ID
        # TODO: the pair ID  (GENE1-GENE2) should be already assigned in the main Ks program
        ids_ = []
        for x in ks.index:
            ids_.append("-".join(sorted([ks.loc[x]['Paralog1'], ks.loc[x]['Paralog2']])))

        ks.index = ids_
        ks_anchors = ks.ix[anchors['pair_id']]

        if out_file:
            ks_anchors.to_csv(out_file)

        return ks, ks_anchors

    else:
        if out_file:
            anchors.to_csv(out_file)

        return anchors


# VISUALIZATION

def stacked_histogram(ks_dist, anchors, ks_range=(0.1, 5), title='$K_S$ distribution', out_file=None):
    """
    Generate a stacked histogram plot.

    :param ks_dist: Ks analysis data frame
    :param anchors: Ks analysis data frame (only anchors)
    :param ks_range:
    :param title:
    :param out_file:
    :return:
    """
    X = ks_dist[ks_dist['Ks'] >= ks_range[0]]
    X = X[X['Ks'] < ks_range[1]]
    A = anchors[anchors['Ks'] >= ks_range[0]]
    A = A[A['Ks'] < ks_range[1]]

    # plot
    plt.figure(figsize=(12, 8))
    plt.title(title, fontsize=18)
    plt.xlabel("Binned $K_S$", fontsize=14)
    plt.hist(X['Ks'], bins=70,
             weights=X['WeightOutliersIncluded'], histtype='barstacked',
             color='black', alpha=0.2, rwidth=0.8, label='Whole paranome')
    plt.hist(A['Ks'], bins=70,
             weights=A['WeightOutliersIncluded'], histtype='barstacked',
             color='black', alpha=0.7, rwidth=0.8, label='Anchors')
    plt.legend(fontsize=12)
    plt.xlim(0.1, 5)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)

    if out_file:
        plt.savefig(out_file, dpi=400)
    else:
        plt.show()


def segments_to_chords_table(segments_file, genome, output_file='chords.tsv'):
    """
    Create chords table for visualization in a chord diagram. Uses the segments.txt
    output of I-ADHoRe. Chords are defined by a source chromosome and a target 
    chromosome and begin and end coordinates for each chromosome respectively.

    TODO: the length of each syntenic block should be included in the table as well
    with length defined as number of genes, not physical length.

    :param segments_file: pat to the I-ADHoRe segments.txt output file
    :param genome: a :func:`gff_parser.Genome object`
    :param output_file: output file name
    """
    segments = pd.read_csv(segments_file, sep='\t', index_col=0)
    segments = segments.groupby('multiplicon')
    chromosomes = segments['list'].apply(list)
    first = segments['first'].apply(list)
    last = segments['last'].apply(list)

    chords = []
    for i in range(1, len(chromosomes) + 1):
        for j in range(len(chromosomes[i])):
            for k in range(j + 1, len(chromosomes[i])):
                source = str(chromosomes[i][j])
                target = str(chromosomes[i][k])

                if source not in genome.genome.keys():
                    logging.warning('Chromosome ID `{}` not found'.format(source))
                elif target not in genome.genome.keys():
                    logging.warning('Chromosome ID `{}` not found'.format(target))
                else:
                    d = {'source_id': source,
                         'source_1': genome.genome[source][first[i][j]]['start'],
                         'source_2': genome.genome[source][last[i][j]]['stop'],
                         'target_id': target,
                         'target_1': genome.genome[target][first[i][k]]['start'],
                         'target_2': genome.genome[target][last[i][k]]['stop'],
                         'label': '{0}-{1};{2}-{3}'.format(first[i][j], last[i][j], first[i][k], last[i][k]),
                         'color': genome.colors[source],
                         'source_length': abs(int(genome.genome[source][last[i][j]]['stop']) -
                                              int(genome.genome[source][first[i][j]]['start']))}
                    chords.append(d)

    df = pd.DataFrame.from_records(chords)
    df.to_csv(output_file, sep='\t')


def visualize(output_dir):
    """
    Visualize intragenomic collinearity
    """
    with open(os.path.join(output_dir, 'vis.html'), 'w') as f:
        f.write(wgd_adhore_html)

    with open(os.path.join(output_dir, 'vis.js'), 'w') as o:
        o.write(circos_js)


# RUNNING EXTERNAL SOFTWARE

def run_adhore(config_file):
    """
    Run I-ADHoRe for a given config file

    :param config_file: path to I-ADHoRe configuration file
    """
    completed = subprocess.run(
        ['i-adhore', config_file], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    logging.warning(completed.stderr.decode('utf-8'))
    logging.info(completed.stdout.decode('utf-8'))
    return


# HTML/JAVSCRIPT TEMPLATES

wgd_adhore_html = """
<!DOCTYPE html>
<meta charset="utf-8">
<head>
    <script src='https://cdn.rawgit.com/nicgirault/circosJS/v2/dist/circos.js'></script>
    <script src="https://d3js.org/d3.v4.min.js"></script>
    <link rel="stylesheet" href="http://www.w3schools.com/lib/w3.css"> 
    <style>
        body {
          font: 10px sans-serif;
        } 
        .ticks {
          font: 10px sans-serif;
        }
        .track,
        .track-inset,
        .track-overlay {
          stroke-linecap: round;
        }
        .track {
          stroke: #000;
          stroke-opacity: 0.3;
          stroke-width: 10px;
        }
        .track-inset {
          stroke: #ddd;
          stroke-width: 8px;
        }
        .track-overlay {
          pointer-events: stroke;
          stroke-width: 50px;
          stroke: transparent;
          cursor: crosshair;
        }
        .handle {
          fill: #fff;
          stroke: #000;
          stroke-opacity: 0.5;
          stroke-width: 1.25px;
        }
    </style>
</head>
<body>

<div class="w3-card-4" style='margin-left:10%;margin-right:40%;margin-bottom:16px;margin-top:16px;'>
	<header class="w3-container w3-green">
  		<h1>Intragenomic collinearity</h1>
	</header>

	<div class="w3-container">
        <p>
		Minimum length (bp) <input style="width:100px;" type="range" min="0" max="1000000" step="1000" value="0">
		</p>
  		<svg id='chart' width='100%', height='800px'></svg>
	</div>

	<footer class="w3-container w3-green">
  		<h5><code>wgd adhore</code> (Arthur Zwaenepoel - 2017)</h5>
	</footer>
</div> 

<div class="w3-card-4" style='margin-left:10%;margin-right:40%;margin-bottom:16px;'>                                                        
    <header class="w3-container w3-green">                                                                               
        <h1><i>K<sub>S</sub></i> distribution</h1>                                                                              
    </header>                                                                                                           
                                                                                                                        
    <div class="w3-container w3-margin">                                                                                          
        <img src='histogram.png' width='100%'>
	</div>                                                                                                              
                                                                                                                        
    <footer class="w3-container w3-green">                                                                               
        <h5><code>wgd adhore</code> (Arthur Zwaenepoel - 2017)</h5>                                                       
    </footer>                                                                                                           
</div>
	
        
    <script>
      var circos = new Circos({
        container: '#chart'
      });
    </script>
    <script src='vis.js'></script>
</body>
</html>
"""

circos_js = """
var minLength = 0;

var drawCircos = function (error, genome, data) {
    var width = 700;
    var circos = new Circos({
        container: '#chart',
        width: width,
        height: width
    })
    
    data = data.filter(function (d) {return parseInt(d.source_length) > minLength})

    data = data.map(function (d) {
        // I think here an if statement can be included for filtering a user defined 
        // syntenic block length
            return {
                source: {
                    id: d.source_id,
                    start: parseInt(d.source_1),
                    end: parseInt(d.source_2),
                    color: d.color,
                    label: d.label
                },
                target: {
                    id: d.target_id,
                    start: parseInt(d.target_1),
                    end: parseInt(d.target_2),
                    color: d.color
                }
            }
        
    })

    circos
        .layout(
            genome,
            {
                innerRadius: width/2 - 80,
                outerRadius: width/2 - 40,
                labels: {
                    radialOffset: 70
                },
                ticks: {
                    display: true,
                    labelDenominator: 1000000
                }
            }
        )
        .chords(
            'l1',
            data,
            {
                opacity: 0.7,
                color: function (d) {return d.source.color;},
                tooltipContent: function (d) {
                    return '<h3>' + d.source.id + ' ➤ ' + d.target.id + ': ' + d.source.label + '</h3><i>(CTRL+C to copy to clipboard)</i>';
                }
            }
        )
        .render()
    }

var svg = d3.select('svg');

d3.queue()
    .defer(d3.json, "genome.json")
    .defer(d3.tsv, "chords.tsv")
    .await(drawCircos);

d3.select("input[type=range]")
    .on("input", inputted);

function inputted() {
      minLength = parseInt(this.value);
      console.log(minLength);
      svg.selectAll("*").remove();
      d3.queue()
        .defer(d3.json, "genome.json")
        .defer(d3.tsv, "chords.tsv")
        .await(drawCircos);
}
"""
