"""
Plotting utilities for wgd
Arthur Zwaenepoel - 2017
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os


def plot_selection(dists, output_file=None, alphas=None, ks_range=(0.1, 5), offset=5, title='Species genus', **kwargs):
    """
    Plot a panel from the `all.csv` output from the Ks analysis.

    :param dists: distribution(s), if multiple provide as list
    :param alphas: alpha values for the different distributions (will assign automatically if not provided)
    :param ks_range: Ks range to include for plotting
    :param offset: offset of axis
    :param title: panel title
    :param kwargs: keyword arguments for :py:func:`matplotlib.pyplot.hist`
    :return: :py:class:`matplotlib.pyplot.Figure` object
    """
    fig = plt.figure(figsize=(12, 12))

    if type(dists) != list:
        dists = [dists]
        alphas = [alphas]

    if not alphas or not alphas[0]:
        alphas = list(np.linspace(0.2, 1, len(dists)))

    for i in range(len(dists)):
        dists[i] = dists[i][dists[i]['Ks'] > ks_range[0]]
        dists[i] = dists[i][dists[i]['Ks'] < ks_range[1]]

    # ks
    ax = fig.add_subplot(221)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([dist['Ks'] for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(dists[i]['Ks'], bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$K_S$')

    # ka
    ax = fig.add_subplot(222)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([dist['Ka'] for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(dists[i]['Ka'], bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$K_A$')

    # log(ka)
    ax = fig.add_subplot(223)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log(dist['Ka']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log(dists[i]['Ka']), bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$ln(K_A)$')

    # log(w)
    ax = fig.add_subplot(224)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log(dist['Omega']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log(dists[i]['Omega']), bins, alpha=alphas[i], color='black', rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$ln(\omega)$')
    fig.suptitle(title)

    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches='tight')

    return fig


def syntenic_dotplot(df, output_file=None):
    """
    Syntenic dotplot function

    :param df: multiplicons pandas data frame
    :param output_file: output file name
    :return: figure
    """
    genomic_elements = {x: 0 for x in list(set(df['list_x']) | set(df['list_y'])) if type(x) == str}

    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111)

    previous = 0
    for key in sorted(genomic_elements.keys()):
        length = max(list(df[df['list_x'] == key]['end_x']) + list(df[df['list_y'] == key]['end_y']))
        genomic_elements[key] = previous
        previous += length

    x = [genomic_elements[key] for key in sorted(genomic_elements.keys())] + [previous]
    ax.vlines(ymin=0, ymax=previous, x=x, linestyles='dotted', alpha=0.2)
    ax.hlines(xmin=0, xmax=previous, y=x, linestyles='dotted', alpha=0.2)
    ax.plot(x, x, color='k', alpha=0.2)
    ax.set_xticks(x)
    ax.set_yticks(x)
    ax.set_xticklabels(x)
    ax.set_yticklabels(x)

    for i in range(len(df)):
        row = df.iloc[i]
        list_x, list_y = row['list_x'], row['list_y']
        if type(list_x) != float:
            curr_list_x = list_x
        x = [genomic_elements[curr_list_x]+x for x in [row['begin_x'], row['end_x']]]
        y = [genomic_elements[list_y]+x for x in [row['begin_y'], row['end_y']]]
        ax.plot(x, y, color='k', alpha=0.5)
        ax.plot(y,x, color='k', alpha=0.5)

    sns.despine(offset=5, trim=True)

    if output_file:
        fig.savefig(output_file, dpi=200, bbox_inches='tight')
        plt.close()

    else:
        return fig


def syntenic_dotplot_ks_colored(df, an, ks, color_map='binary', output_file=None):
    """
    Syntenic dotplot with segment colored by mean Ks value

    :param df: multiplicons pandas data frame
    :param an: anchorpoints pandas data frame
    :param ks: Ks distribution data frame
    :param color_map: color map string
    :param output_file: output file name
    :return: figure
    """
    cmap = plt.get_cmap(color_map)

    an['pair'] = an['gene_x'].astype(str) + '-' + an['gene_y']
    an['pair'] = an['pair'].map(lambda x: '-'.join(sorted(x.split('-'))))
    genomic_elements = {x: 0 for x in list(set(df['list_x']) | set(df['list_y'])) if type(x) == str}
    ks_multiplicons = {}
    all_ks = []
    for i in range(len(df)):
        row = df.iloc[i]
        pairs = an[an['multiplicon'] == row['id']]['pair']
        mean_ks = np.mean(ks.loc[pairs]['Ks'])
        ks_multiplicons[row['id']] = mean_ks
        if mean_ks < 5:
            all_ks.append(mean_ks)
    z = [[0, 0], [0, 0]]
    levels = range(0, 101, 1)
    tmp = plt.contourf(z, levels, cmap=cmap)
    plt.clf()

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)

    previous = 0
    for key in sorted(genomic_elements.keys()):
        length = max(list(df[df['list_x'] == key]['end_x']) + list(df[df['list_y'] == key]['end_y']))
        genomic_elements[key] = previous
        previous += length

    x = [genomic_elements[key] for key in sorted(genomic_elements.keys())] + [previous]
    ax.vlines(ymin=0, ymax=previous, x=x, linestyles='dotted', alpha=0.2)
    ax.hlines(xmin=0, xmax=previous, y=x, linestyles='dotted', alpha=0.2)
    ax.plot(x, x, color='k', alpha=0.2)
    ax.set_xticks(x)
    ax.set_yticks(x)
    ax.set_xticklabels(x)
    ax.set_yticklabels(x)

    for i in range(len(df)):
        row = df.iloc[i]
        list_x, list_y = row['list_x'], row['list_y']
        if type(list_x) != float:
            curr_list_x = list_x
        x = [genomic_elements[curr_list_x] + x for x in [row['begin_x'], row['end_x']]]
        y = [genomic_elements[list_y] + x for x in [row['begin_y'], row['end_y']]]
        ax.plot(x, y, alpha=0.7, linewidth=3, color=cmap(ks_multiplicons[row['id']] / 5))
        ax.plot(y, x, alpha=0.7, linewidth=3, color=cmap(ks_multiplicons[row['id']] / 5))

    cbar = plt.colorbar(tmp, fraction=0.02, pad=0.01)
    cbar.ax.set_yticklabels(['{:.2f}'.format(x) for x in np.linspace(0, 5, 11)])
    sns.despine(offset=5, trim=True)

    if output_file:
        fig.savefig(output_file, dpi=200, bbox_inches='tight')
        plt.close()

    else:
        return fig


def visualize_adhore_circos_js(output_dir):
    """
    Visualize intragenomic collinearity
    """
    with open(os.path.join(output_dir, 'vis.html'), 'w') as f:
        f.write(wgd_adhore_html)

    with open(os.path.join(output_dir, 'vis.js'), 'w') as o:
        o.write(circos_js)


# HTML/JAVSCRIPT TEMPLATES ---------------------------------------------------------------------------------------------

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
