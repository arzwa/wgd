"""
Arthur Zwaenepoel - 2017
"""
import plumbum as pb
import matplotlib
if not 'DISPLAY' in pb.local.env:
    matplotlib.use('Agg')  # use this backend when no X server
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
import os
import matplotlib.patheffects as pe
import pandas as pd
import random
from .modeling import weighted_to_unweighted


def plot_selection(dists, output_file=None, alphas=None, colors=None, labels=None, ks_range=(0.1, 5), offset=5,
                   title='Species genus', **kwargs):
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
    fig = plt.figure(figsize=(12, 10))

    if type(dists) != list:
        dists = [dists]
        alphas = [alphas]

    if not alphas or not alphas[0]:
        alphas = list(np.linspace(0.2, 1, len(dists)))

    if not colors or not colors[0]:
        colors = ['black']*len(dists)

    if not labels or not labels[0]:
        labels = [None] * len(dists)

    for i in range(len(dists)):
        dists[i] = dists[i][dists[i]['Ks'] > ks_range[0]]
        dists[i] = dists[i][dists[i]['Ks'] < ks_range[1]]

    # ks
    ax = fig.add_subplot(221)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([dist['Ks'] for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(dists[i]['Ks'], bins, alpha=alphas[i], color=colors[i], rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], label=labels[i], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$K_S$')

    # ka
    ax = fig.add_subplot(222)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log10(dist['Ks']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log10(dists[i]['Ks']), bins, alpha=alphas[i], color=colors[i], rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], label=labels[i], **kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$log_{10}(K_s)$')

    # log(ka)
    ax = fig.add_subplot(223)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log10(dist['Ka']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log10(dists[i]['Ka']), bins, alpha=alphas[i], color=colors[i], rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], label=labels[i],**kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$log_{10}(K_A)$')

    # log(w)
    ax = fig.add_subplot(224)
    # get the bin edges
    bins = np.histogram(np.hstack(tuple([np.log10(dist['Omega']) for dist in dists])), bins=40)[1]
    for i in range(len(dists)):
        ax.hist(np.log10(dists[i]['Omega']), bins, alpha=alphas[i], color=colors[i], rwidth=0.8,
                weights=dists[i]['WeightOutliersIncluded'], label=labels[i],**kwargs)
    sns.despine(offset=offset, trim=True)
    ax.set_xlabel('$log_{10}(\omega)$')

    if labels[0]:
        plt.legend()

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
    genomic_elements_ = {x: 0 for x in list(set(df['list_x']) | set(df['list_y'])) if type(x) == str}

    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111)

    for key in sorted(genomic_elements_.keys()):
        length = max(list(df[df['list_x'] == key]['end_x']) + list(df[df['list_y'] == key]['end_y']))
        genomic_elements_[key] = length

    previous = 0
    genomic_elements = {}
    sorted_ge = sorted(genomic_elements_.items(), key=lambda x: x[1], reverse=True)
    for kv in sorted_ge:
        genomic_elements[kv[0]] = previous
        previous += kv[1]

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


def syntenic_dotplot_ks_colored(df, an, ks, color_map='Spectral', output_file=None):
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
    genomic_elements_ = {x: 0 for x in list(set(df['list_x']) | set(df['list_y'])) if type(x) == str}
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

    fig = plt.figure(figsize=(16,15))
    ax = fig.add_subplot(111)

    for key in sorted(genomic_elements_.keys()):
        length = max(list(df[df['list_x'] == key]['end_x']) + list(df[df['list_y'] == key]['end_y']))
        genomic_elements_[key] = length

    previous = 0
    genomic_elements = {}
    sorted_ge = sorted(genomic_elements_.items(), key=lambda x: x[1], reverse=True)
    for kv in sorted_ge:
        genomic_elements[kv[0]] = previous
        previous += kv[1]

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
        ax.plot(x, y, alpha=0.7, linewidth=3, color=cmap(ks_multiplicons[row['id']] / 5),
                path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
        ax.plot(y, x, alpha=0.7, linewidth=3, color=cmap(ks_multiplicons[row['id']] / 5),
                path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

    cbar = plt.colorbar(tmp, fraction=0.02, pad=0.01)
    cbar.ax.set_yticklabels(['{:.2f}'.format(x) for x in np.linspace(0, 5, 11)])
    sns.despine(offset=5, trim=True)

    if output_file:
        fig.savefig(output_file, dpi=200, bbox_inches='tight')
        plt.close()

    else:
        return fig


def histogram_bokeh(ks_distributions, labels):
    # TODO
    from bokeh.io import curdoc
    from bokeh.layouts import row, column
    from bokeh.models import ColumnDataSource
    from bokeh.models.widgets import PreText, Select, TextInput, MultiSelect, Slider
    from bokeh.plotting import figure
    from bokeh.client import push_session
    from pylab import cm, colors
    from scipy.stats import gaussian_kde as kde

    def get_colors(cmap_choice='hsv'):
        cmap = cm.get_cmap(cmap_choice, len(ks_distributions))
        if cmap_choice == 'binary':
            cmap = cm.get_cmap(cmap_choice, len(ks_distributions)*2)
        c = []
        for i in range(cmap.N):
            rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
            c.append(matplotlib.colors.rgb2hex(rgb))
        if cmap_choice == 'binary':
            c = c[len(ks_distributions):]
        return c

    def get_data(df, var, scale, r1, r2):
        df = df[df[var] > r1]
        df = df[df[var] < r2]
        data = df[var].dropna()
        if scale == 'log10':
            data = np.log10(data)
        weights = df['WeightOutliersIncluded']
        return data, weights, df

    dists = [pd.read_csv(x, sep='\t') for x in ks_distributions]
    if labels:
        labels = labels.split(',')
    else:
        labels= ks_distributions
    dist_ids = {ks_distributions[i]: i for i in range(len(ks_distributions))}
    df = dists[0]
    c = get_colors()
    print(c)
    variables = ['Ks', 'Ka', 'Omega']
    scales = ['Normal', 'log10']

    stats = PreText(text='', width=500)
    var = Select(title='Variable', value='Ks', options=variables)
    scale = Select(title='Scale', value='Normal', options=scales)
    r1 = TextInput(title="Minimum", value='0.1')
    r2 = TextInput(title="Maximum", value='5')
    bins = TextInput(title="Bins", value='50')
    line = Slider(title="Lines", start=0, end=1, value=0.3, step=0.1)
    density = Slider(title="Kernel density", start=0, end=1, value=0, step=1)
    selected_dists = MultiSelect(title='Distributions', options=ks_distributions, value=[ks_distributions[0]])
    color_choice = Select(options=['binary', 'hsv', 'hot', 'magma', 'viridis', 'Greens'],
                          value='hsv', title='Color map')

    data, weights, df = get_data(df, var.value, scale.value, float(r1.value), float(r2.value))
    all_data = [data]

    p1 = figure(plot_width=1000, plot_height=700, tools='pan,wheel_zoom,box_select,reset')
    p1.xgrid.grid_line_color = None
    p1.ygrid.grid_line_color = None
    p1.border_fill_color = 'white'
    p1.outline_line_color = None
    p1.yaxis.axis_label = '# paralogs'
    p1.xaxis.axis_label = 'Ks'

    hist, edges = np.histogram(data, bins=int(bins.value), weights=weights)
    hist_dict = {0: p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color=c[0],
                            line_color=c[0], fill_alpha=0.4, line_alpha=line.value, legend=labels[0])}
    p1.legend.label_text_font_style = "italic"
    density_dict = {}
    # set up callbacks

    def update(selected=None):
        all_data = []
        c = get_colors(color_choice.value)
        for f in selected_dists.value:
            i = dist_ids[f]
            df = dists[i]
            data, weights, df = get_data(df, var.value, scale.value, float(r1.value), float(r2.value))
            hist, edges = np.histogram(data, bins=int(bins.value), weights=weights)
            if i in hist_dict:
                remove_plot(hist_dict, i)
            hist_dict[i] = p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color=c[i],
                                       line_color=c[i], fill_alpha=0.4, line_alpha=line.value,
                                       legend=labels[i])

            if density.value == 1:
                dw = [x[0] for x in weighted_to_unweighted(df)]
                data_weighed = np.random.choice(dw, 10000)
                dens= kde(data_weighed)
                xs = np.arange(float(r1.value), float(r2.value), 0.01)
                ys = dens(xs)
                conv_factor = min(hist)/min(ys)
                ys *= conv_factor
                density_dict[i] = p1.line(xs, ys, color=c[i])

            x_label = var.value
            if scale.value == 'log10':
                x_label = 'log10(' + var.value + ')'
            p1.xaxis.axis_label = x_label

            all_data.append(data)

        for d in dist_ids:
            i = dist_ids[d]
            if d not in selected_dists.value and i in hist_dict.keys():
                print('updating {}'.format(i))
                remove_plot(hist_dict, i)
                #hist_dict[i].data_source.data['line_alpha'] = []

        update_stats(all_data)

    def update_stats(data):
        df = pd.concat(data, axis=1)
        stats.text = str(df.describe())

    def remove_plot(hist_dict, i):
        hist_dict[i].data_source.data["top"] = []
        hist_dict[i].data_source.data["left"] = []
        hist_dict[i].data_source.data["right"] = []

    def change_update(attrname, old, new):
        update()

    var.on_change('value', change_update)
    scale.on_change('value', change_update)
    r1.on_change('value', change_update)
    r2.on_change('value', change_update)
    line.on_change('value', change_update)
    bins.on_change('value', change_update)
    color_choice.on_change('value', change_update)
    density.on_change('value', change_update)
    selected_dists.on_change('value', change_update)

    # set up layout
    widgets1 = column(selected_dists, var, scale, r1, r2, bins, line, density, color_choice, stats)
    main_row = row(p1, widgets1)
    layout = column(main_row)

    # initialize
    session = push_session(curdoc())
    #update()
    curdoc().add_root(layout)
    #curdoc().add_periodic_callback(update, 50)
    session.show(layout)  # open the document in a browser
    session.loop_until_closed()  # run forever
    pass


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
                    return '<h3>' + d.source.id + ' > ' + d.target.id + ': ' + d.source.label + '</h3><i>(CTRL+C to copy to clipboard)</i>';
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
