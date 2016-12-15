import numpy as np
import matplotlib
matplotlib.use("Agg")
try:
    from scipy.stats.kde import gaussian_kde
    CAN_PLOT_KDE = True
except ImportError:
    print "KDE will be omitted from window distribution plot due to missing scipy library."
    CAN_PLOT_KDE = False

import pylab
matplotlib.rc("mathtext", fontset="stix")

LINE_WIDTH = 1.25 # table borders and ticks
MARKER_SIZE = 5

FONT_PARAMS = dict(
    family='serif',
    serif='Times New Roman',
    size=12,
    )
x_label_size_increase = 4
y_label_size_increase = 1

matplotlib.rc('font', **FONT_PARAMS)
matplotlib.rc('axes', linewidth=LINE_WIDTH)

def create_figure(figsize=(5,4)):
    fig = matplotlib.pyplot.figure(figsize=figsize)
    fig.hold = True
    return fig

def add_axis_to_figure(fig, subplot_layout=111, sharex=None, sharey=None):
    return fig.add_subplot(subplot_layout, sharex=sharex, sharey=sharey)

def save_figure(fig, fig_name):
    fig.savefig(fig_name)

def simple_plot(xs, ys, es, xlabel, ylabel, fig_name=None, show_auc=False):
    fig = create_figure()
    ax = add_axis_to_figure(fig)
    plot(ax, xs, ys, es, xlabel=xlabel, ylabel=ylabel, color="b", show_auc=show_auc)

    finalise_figure(fig, ax, xlabel, ylabel, fig_name)

def plot(ax, xs, ys, es, xlabel=None, ylabel=None, label=None, color="b", linewidth=LINE_WIDTH, alpha=1, line_style="-", show_auc=False):

    ax.errorbar(xs, ys, yerr=es, marker="o", zorder=1, linestyle=line_style, color=color, markerfacecolor=color, label=label, markersize=MARKER_SIZE, alpha=alpha, linewidth=linewidth)
    if show_auc:
        ax.fill_between(xs, ys, np.zeros(len(xs)), facecolor='blue', alpha=0.2)

def plot_with_kde(ax, xs, ys, xlabel=None, ylabel=None, label=None, show=True, kde=None, color="blue"):
    ax.plot(xs, ys, color=color, label=label)
    if kde:
        kde_curve = kde(np.array(xs))
        ax.fill_between(xs, kde_curve, alpha=0.5, facecolor=color)

def plot_hist_to_axis(ax, centers, his, kde=None):
    plot_with_kde(ax, centers, his, show=False, kde=kde)

def finalise_figure(fig, ax, xlabel, ylabel, fig_name=None, show=False):
    pylab.ylabel(ylabel, fontsize=FONT_PARAMS["size"] + y_label_size_increase)
    pylab.xlabel(xlabel, fontsize=FONT_PARAMS["size"] + x_label_size_increase)
    pylab.xticks(fontsize=FONT_PARAMS["size"])
    pylab.yticks(fontsize=FONT_PARAMS["size"])

    [i.set_linewidth(LINE_WIDTH) for i in ax.spines.itervalues()]
    fig.tight_layout()
    if fig_name:
        save_figure(fig, fig_name)
    if show:
        pylab.show()
    matplotlib.pyplot.close(fig)

def generate_histogram(values, n_bins=50):
    his, bins = np.histogram(values, bins=n_bins, normed=True)
    centers = (bins[:-1]+bins[1:])/2
    kde = gaussian_kde(values) if CAN_PLOT_KDE else None
    return centers, his, kde
