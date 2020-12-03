# TODO: still needs to be updated
import plumbum as pb
import matplotlib

if not 'DISPLAY' in pb.local.env:
    matplotlib.use('Agg')  # use this backend when no X server
import matplotlib.pyplot as plt
import matplotlib
import logging
import numpy as np
import seaborn as sns
import pandas as pd


def node_averages(df):
    # note that this returns a df with fewer rows, i.e. one for every
    # node in the gene family trees.
    return df.groupby(["family", "node"])["dS"].mean()

def node_weights(df):
    # note that this returns a df with the same number of rows
    return 1 / df.groupby(["family", "node"])["dS"].transform('count')

def parse_filter(s):
    x = [x.strip() for x in s.split("<")]
    if len(x) != 3:
        raise(ValueError("invalid 'x < field < y' filter string"))
    return (x[1], float(x[0]), float(x[2]))

def parse_filters(filterstring):
    return [parse_filter(s) for s in filterstring.split(",")]

def apply_filters(df, filters):
    for key, lower, upper in filters:
        df = df[df[key] > lower]
        df = df[df[key] < upper]
    return df

_labels = {
        "dS" : "$K_\mathrm{S}$",
        "dN" : "$K_\mathrm{A}$",
        "dN/dS": "$\omega$"}

def default_plot(
        *args, 
        alphas=None,
        colors=None,
        weighted=True, 
        title="",
        ylabel="duplication events",
        **kwargs):
    """
    Make a figure of node-weighted histograms for multiple distributions and
    variables. Returns the figure object.
    
    !!! note: Assumes the data frames are filtered as desired. 
    """
    ndists = len(args)
    alphas = alphas or list(np.linspace(0.2, 1, ndists))
    colors = colors or ['black'] * ndists
    
    # assemble panels
    keys = ["dS", "dS", "dN", "dN/dS"]
    funs = [lambda x: x, np.log10, np.log10, np.log10]
    fig, axs = plt.subplots(2, 2)
    for (c, a, dist) in zip(colors, alphas, args):
        for ax, k, f in zip(axs.flatten(), keys, funs):
            w = node_weights(dist)
            x = f(dist[k])
            y = x[np.isfinite(x)]
            w = w[np.isfinite(x)]
            ax.hist(y, weights=w, color=c, alpha=a, **kwargs)
            xlabel = _labels[k]
            if f == np.log10:
                xlabel = "$\log_{10}" + xlabel[1:-1] + "$"
            ax.set_xlabel(xlabel)
    axs[0,0].set_ylabel(ylabel)
    axs[1,0].set_ylabel(ylabel)

    # finalize plot
    sns.despine(offset=1)
    fig.suptitle(title, x=0.1, y=0.99, ha="left", va="top")
    fig.tight_layout()
    plt.subplots_adjust(top=0.85)  # prevent suptitle from overlapping
    return fig


#def syntenic_dotplot(df, min_length=250, output_file=None):
#    """
#    Syntenic dotplot function
#
#    :param df: multiplicons pandas data frame
#    :param min_length: minimum length of a genomic element
#    :param output_file: output file name
#    :return: figure
#    """
#    genomic_elements_ = {
#        x: 0 for x in list(set(df['list_x']) | set(df['list_y']))
#        if type(x) == str
#    }
#
#    fig = plt.figure(figsize=(6.5, 6))
#    ax = fig.add_subplot(111)
#
#    for key in sorted(genomic_elements_.keys()):
#        length = max(list(df[df['list_x'] == key]['end_x']) + list(
#                df[df['list_y'] == key]['end_y']))
#        if length >= min_length:
#            genomic_elements_[key] = length
#
#    previous = 0
#    genomic_elements = {}
#    sorted_ge = sorted(genomic_elements_.items(), key=lambda x: x[1],
#                       reverse=True)
#    labels = [kv[0] for kv in sorted_ge if kv[1] >= min_length]
#
#    for kv in sorted_ge:
#        genomic_elements[kv[0]] = previous
#        previous += kv[1]
#
#    x = [genomic_elements[key] for key in sorted(genomic_elements.keys())] + \
#        [previous]
#    x = sorted(list(set(x)))  # FIXME hack
#    if len(x) == 0:
#        logging.warning("No multiplicons found!")
#        return
#
#    # plot layout stuff!
#    ax.vlines(ymin=0, ymax=previous, x=x, linestyles='dotted', alpha=0.2)
#    ax.hlines(xmin=0, xmax=previous, y=x, linestyles='dotted', alpha=0.2)
#    ax.plot(x, x, color='k', alpha=0.2)
#    ax.set_xticks(x)
#    ax.set_yticks(x)
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
#    ax.set_xlim(0, max(x))
#    ax.set_ylim(0, max(x))
#    ax.set_xticks([(x[i] + x[i - 1]) / 2 for i in range(1, len(x))], minor=True)
#    ax.set_xticklabels(labels, minor=True, rotation=45)
#    ax.set_yticks([(x[i] + x[i - 1]) / 2 for i in range(1, len(x))], minor=True)
#    ax.set_yticklabels(labels, minor=True, rotation=45)
#    for tick in ax.xaxis.get_minor_ticks():
#        tick.tick1line.set_markersize(0)
#        tick.tick2line.set_markersize(0)
#        tick.label1.set_horizontalalignment('center')
#    for tick in ax.yaxis.get_minor_ticks():
#        tick.tick1line.set_markersize(0)
#        tick.tick2line.set_markersize(0)
#
#    # actual dots (or better, line segments)
#    for i in range(len(df)):
#        row = df.iloc[i]
#        list_x, list_y = row['list_x'], row['list_y']
#        if type(list_x) != float:
#            curr_list_x = list_x
#        x = [genomic_elements[curr_list_x] + x for x in
#             [row['begin_x'], row['end_x']]]
#        y = [genomic_elements[list_y] + x for x in
#             [row['begin_y'], row['end_y']]]
#        ax.plot(x, y, color='k', alpha=0.7)
#        ax.plot(y, x, color='k', alpha=0.7)
#
#    # saving
#    if output_file:
#        fig.savefig(output_file, dpi=200, bbox_inches='tight')
#        plt.close()
#
#    else:
#        return fig
#
#
#def syntenic_dotplot_ks_colored(
#        df, an, ks, min_length=50, color_map='Spectral', min_ks=0.05, max_ks=5,
#        output_file=None
#):
#    """
#    Syntenic dotplot with segment colored by mean Ks value
#
#    :param df: multiplicons pandas data frame
#    :param an: anchorpoints pandas data frame
#    :param ks: Ks distribution data frame
#    :param min_length: minimum length of a genomic element
#    :param color_map: color map string
#    :param min_ks: minimum median Ks value
#    :param max_ks: maximum median Ks value
#    :param output_file: output file name
#    :return: figure
#    """
#    cmap = plt.get_cmap(color_map)
#    if len(an["gene_x"]) == 0:
#        logging.warning("No multiplicons found!")
#        return
#    an["pair"] = an.apply(lambda x: '__'.join(
#            sorted([x["gene_x"], x["gene_y"]])), axis=1)
#    genomic_elements_ = {
#        x: 0 for x in list(set(df['list_x']) | set(df['list_y']))
#        if type(x) == str
#    }
#
#    ks_multiplicons = {}
#    all_ks = []
#    for i in range(len(df)):
#        row = df.iloc[i]
#        pairs = an[an['multiplicon'] == row['id']]['pair']
#        med_ks = np.median(ks.loc[ks.index.intersection(pairs)]['Ks'])
#        ks_multiplicons[row['id']] = med_ks
#        all_ks.append(med_ks)
#
#    z = [[0, 0], [0, 0]]
#    levels = range(0, 101, 1)
#    tmp = plt.contourf(z, levels, cmap=cmap)
#    plt.clf()
#
#    fig = plt.figure(figsize=(6.5, 6))
#    ax = fig.add_subplot(111)
#
#    for key in sorted(genomic_elements_.keys()):
#        length = max(list(df[df['list_x'] == key]['end_x']) + list(
#                df[df['list_y'] == key]['end_y']))
#        if length >= min_length:
#            genomic_elements_[key] = length
#
#    previous = 0
#    genomic_elements = {}
#    sorted_ge = sorted(genomic_elements_.items(), key=lambda x: x[1],
#                       reverse=True)
#    labels = [kv[0] for kv in sorted_ge if kv[1] >= min_length]
#
#    for kv in sorted_ge:
#        genomic_elements[kv[0]] = previous
#        previous += kv[1]
#
#    # plot layout
#    x = [genomic_elements[key] for key in sorted(genomic_elements.keys())] + \
#        [previous]
#    x = sorted(list(set(x)))
#    ax.vlines(ymin=0, ymax=previous, x=x, linestyles='dotted', alpha=0.2)
#    ax.hlines(xmin=0, xmax=previous, y=x, linestyles='dotted', alpha=0.2)
#    ax.plot(x, x, color='k', alpha=0.2)
#    ax.set_xticks(x)
#    ax.set_yticks(x)
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
#    ax.set_xlim(0, max(x))
#    ax.set_ylim(0, max(x))
#    ax.set_xticks([(x[i] + x[i - 1]) / 2 for i in range(1, len(x))], minor=True)
#    ax.set_xticklabels(labels, minor=True, rotation=45)
#    ax.set_yticks([(x[i] + x[i - 1]) / 2 for i in range(1, len(x))], minor=True)
#    ax.set_yticklabels(labels, minor=True, rotation=45)
#    for tick in ax.xaxis.get_minor_ticks():
#        tick.tick1line.set_markersize(0)
#        tick.tick2line.set_markersize(0)
#        tick.label1.set_horizontalalignment('center')
#    for tick in ax.yaxis.get_minor_ticks():
#        tick.tick1line.set_markersize(0)
#        tick.tick2line.set_markersize(0)
#        # tick.label1.set_horizontalalignment('center')
#
#    # the actual dots (or better, line segments)
#    for i in range(len(df)):
#        row = df.iloc[i]
#        list_x, list_y = row['list_x'], row['list_y']
#        if type(list_x) != float:
#            curr_list_x = list_x
#        x = [genomic_elements[curr_list_x] + x for x in
#             [row['begin_x'], row['end_x']]]
#        y = [genomic_elements[list_y] + x for x in
#             [row['begin_y'], row['end_y']]]
#        med_ks = ks_multiplicons[row['id']]
#        if min_ks < med_ks <= max_ks:
#            ax.plot(x, y, alpha=0.9, linewidth=3,
#                    color=cmap(ks_multiplicons[row['id']] / 5)),
#            # path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
#            ax.plot(y, x, alpha=0.9, linewidth=3,
#                    color=cmap(ks_multiplicons[row['id']] / 5))
#            # path_effects=[pe.Stroke(linewidth=4, foreground='k'),
#            # pe.Normal()])
#
#    # colorbar
#    cbar = plt.colorbar(tmp, fraction=0.02, pad=0.01)
#    cbar.ax.set_yticklabels(['{:.2f}'.format(x) for x in np.linspace(0, 5, 11)])
#
#    # saving
#    if output_file:
#        fig.savefig(output_file, dpi=200, bbox_inches='tight')
#        plt.close()
#
#    else:
#        return fig
#
#
#def histogram_bokeh(ks_distributions, labels):
#    """
#    Run an interactive bokeh application.
#    This requires a running bokeh server! Use ``bokeh serve &`` to start a bokeh
#    server in the background.
#
#    :param ks_distributions: a list of Ks distributions (pandas data frames)
#    :param labels: a list of labels for the corresponding distributions
#    :return: bokeh app
#    """
#    from bokeh.io import curdoc
#    from bokeh.layouts import widgetbox, layout
#    from bokeh.models.widgets import Select, TextInput, Slider, Div
#    from bokeh.models.widgets import CheckboxGroup, Toggle
#    from bokeh.plotting import figure, output_file, show
#    from bokeh.client import push_session
#    from pylab import cm
#    from .utils import gaussian_kde
#    from .modeling import reflect
#
#    # helper functions
#    def get_colors(cmap_choice='binary'):
#        too_light = [
#            'binary', 'hot', 'copper', 'pink', 'summer', 'bone', 'Pastel1',
#            'Pastel2', 'gist_ncar', 'nipy_spectral', 'Greens'
#        ]
#        cmap = cm.get_cmap(cmap_choice, len(ks_distributions))
#        if cmap_choice in too_light:
#            cmap = cm.get_cmap(cmap_choice, len(ks_distributions) * 4)
#        c = []
#        for i in range(cmap.N):
#            rgb = cmap(i)[:3]  # will return rgba, but we need only rgb
#            c.append(matplotlib.colors.rgb2hex(rgb))
#        if cmap_choice in too_light:
#            if len(ks_distributions) > 1:
#                c = c[2:-2]
#            else:
#                c = [c[-1]]
#        return c
#
#    def get_data(df, var, scale, r1, r2, outliers_included):
#        df = filter_group_data(df, min_ks=r1, max_ks=r2,
#                               weights_outliers_included=outliers_included)
#        data = df[var].dropna()
#        if scale == 'log10':
#            data = np.log10(data)
#        return data, df
#
#    # get the distributions
#    dists = [
#        pd.read_csv(x, sep='\t')
#        for x in ks_distributions
#    ]
#    if labels:
#        labels = labels.split(',')
#    else:
#        labels = ks_distributions
#
#    # basics
#    c = get_colors()
#    variables = ['Ks', 'Ka', 'Omega']
#    scales = ['Normal', 'log10']
#
#    # set up widgets
#    div = Div(text=BOKEH_APP_DIV, width=800)
#    var = Select(title='Variable', value='Ks', options=variables)
#    scale = Select(title='Scale', value='Normal', options=scales)
#    r1 = TextInput(title="Minimum", value='0.1')
#    r2 = TextInput(title="Maximum", value='5')
#    bins = TextInput(title="Bins", value='50')
#    bandwidth = TextInput(title="Bandwidth", value='0.1')
#    line = Slider(title="Lines", start=0, end=1, value=0.3, step=0.1)
#    density = CheckboxGroup(labels=["Histogram", "KDE"], active=[0])
#    density_alpha = Slider(title="Density alpha value", start=0, end=1,
#                           value=0.6, step=0.1)
#    hist_alpha = Slider(title="Histogram alpha value", start=0, end=1,
#                        value=0.6, step=0.1)
#    color_choice = Select(
#            options=[
#                'binary', 'hsv', 'hot', 'magma', 'viridis', 'Greens',
#                'spring', 'autumn', 'copper', 'cool', 'winter', 'pink',
#                'summer', 'bone', 'RdBu', 'RdYlGn', 'coolwarm', 'inferno',
#                'Pastel1', 'Pastel2', 'tab10', 'gnuplot', 'brg', 'gist_ncar',
#                'jet', 'rainbow', 'nipy_spectral', 'ocean', 'cubehelix'
#            ],
#            value='binary', title='Color map'
#    )
#    no_reweight = Toggle(
#            label="Don't adapt weights when filtering", active=False)
#
#    # set up figure
#    p1 = figure(plot_width=1000, plot_height=700,  # output_backend="svg",
#                tools='pan,wheel_zoom,xwheel_zoom,ywheel_zoom,save')
#    p1.xgrid.grid_line_color = None
#    p1.ygrid.grid_line_color = None
#    p1.border_fill_color = 'white'
#    p1.outline_line_color = None
#    p1.yaxis.axis_label = 'Duplications'
#    p1.xaxis.axis_label = 'Ks'
#
#    # draw initial plot
#    hist_dict = {}
#    density_dict = {}
#    all_data = []
#
#    # set up callbacks
#    def update(selected=None):
#        redraw_plots()
#
#    def redraw_plots():
#        print(density.active)
#        c = get_colors(color_choice.value)
#        p1.legend.items = []
#
#        all_data = []
#        for i in range(len(dists)):
#            df = dists[i]
#            data, df = get_data(
#                    df, var.value, scale.value, float(r1.value),
#                    float(r2.value), no_reweight.active
#            )
#            all_data.append(data)
#
#        edges = np.histogram(
#                np.hstack(tuple(all_data)), bins=int(bins.value))[1]
#
#        for i in range(len(dists)):
#            if density.active == [0]:
#                hist = np.histogram(all_data[i], bins=int(bins.value))[0]
#                p1.yaxis.axis_label = 'Duplications'
#            else:
#                hist = np.histogram(
#                        all_data[i], bins=int(bins.value), density=True)[0]
#                p1.yaxis.axis_label = 'density'
#
#            # First histograms
#            if i in hist_dict:
#                remove_plot(hist_dict, i)
#
#            if 0 in density.active:
#                hist_dict[i] = p1.quad(
#                        top=hist, bottom=0, left=edges[:-1],
#                        right=edges[1:], fill_color=c[i],
#                        line_color="black",
#                        fill_alpha=hist_alpha.value,
#                        line_alpha=line.value,
#                        legend=labels[i]
#                )
#
#            # Then KDEs
#            if i in density_dict:
#                density_dict[i].data_source.data['x'] = []
#                density_dict[i].data_source.data['y'] = []
#
#            if 1 in density.active:
#                X = reflect(all_data[i])
#                kde = gaussian_kde(X, bw_method=float(bandwidth.value))
#                x = np.linspace(
#                        float(r1.value) + 0.000001, float(r2.value), 1000)
#                if scale.value == 'log10':
#                    x = np.log10(x)
#                pdf = np.array(kde(x)) * 2
#
#                # add boundaries such that it is a nice curve!
#                pdf = np.hstack([0, pdf, 0])
#                if scale.value == 'log10':
#                    x = np.hstack([
#                        np.log10(float(r1.value) + 0.00000099), x,
#                        np.log10(float(r2.value) + 0.000001)
#                    ])
#                else:
#                    x = np.hstack(
#                            [float(r1.value), x, float(r2.value) + 0.000001])
#
#                density_dict[i] = p1.patch(
#                        x=x, y=pdf, fill_color=c[i], line_width=2,
#                        line_color="black", line_alpha=line.value,
#                        alpha=density_alpha.value, legend=labels[i]
#                )
#
#            p1.legend.label_text_font_style = "italic"
#            p1.legend.click_policy = "hide"
#            p1.legend.inactive_fill_alpha = 0.6
#            v = var.value
#            if v == "Omega":
#                v = "Ka/Ks"
#            if scale.value == 'log10':
#                v = 'log10(' + v + ')'
#            p1.xaxis.axis_label = v
#
#    def remove_plot(hist_dict, i):
#        hist_dict[i].data_source.data["top"] = []
#        hist_dict[i].data_source.data["left"] = []
#        hist_dict[i].data_source.data["right"] = []
#
#    def change_update(attrname, old, new):
#        update()
#
#    def feat_change(attrname, old, new):
#        c = get_colors(color_choice.value)
#        for i, d in density_dict.items():
#            d.glyph.fill_color = c[i]
#            d.glyph.line_color = "black"
#            d.glyph.fill_alpha = density_alpha.value
#            d.glyph.line_alpha = line.value
#        for i, d in hist_dict.items():
#            d.glyph.fill_color = c[i]
#            d.glyph.line_color = "black"
#            d.glyph.fill_alpha = hist_alpha.value
#            d.glyph.line_alpha = line.value
#
#    def bins_update(attrname, old, new):
#        update()
#
#    var.on_change('value', bins_update)
#    bandwidth.on_change('value', bins_update)
#    scale.on_change('value', bins_update)
#    r1.on_change('value', bins_update)
#    r2.on_change('value', bins_update)
#    line.on_change('value', feat_change)
#    bins.on_change('value', bins_update)
#    color_choice.on_change('value', feat_change)
#    density.on_change('active', bins_update)
#    density_alpha.on_change('value', feat_change)
#    hist_alpha.on_change('value', feat_change)
#    no_reweight.on_change("active", bins_update)
#
#    # set up layout
#    widgets1 = widgetbox(
#            var, scale, color_choice, density, line, hist_alpha,
#            density_alpha, r1, r2, bins, bandwidth, no_reweight,
#            sizing_mode='fixed'
#    )
#    l = layout([
#        [div],
#        [widgets1, p1],
#    ], sizing_mode='fixed')
#
#    # initialize
#    update()
#
#    session = push_session(curdoc())
#    curdoc().add_root(l)
#    session.show(l)  # open the document in a browser
#    session.loop_until_closed()  # run forever
#    return  # show(p1)
#
#
## HTML/JAVSCRIPT TEMPLATES -----------------------------------------------------
#BOKEH_APP_DIV = """
#<div>
#    <h1><i>K</i><sub>S</sub> distribution visualization</h1>
#    <p> 
#        Use the widgets on the right to tune the plot in accordance with your 
#        desires. You can click on the legend names to hide/show the relevant 
#        distribution. Note that for more sensible legend names (by default the 
#        file names are used), the <code>--label / -l</code> flag can come in 
#        handy. Note that KDEs are corrected for left boundary effects by 
#        reflection around the minimum <i>K</i><sub>S</sub> value.
#    </p>
#</div>
#"""
