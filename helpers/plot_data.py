import os
import matplotlib
if not os.environ.has_key("DISPLAY"):
    matplotlib.use("Agg")
import pylab
import numpy as np
try:
    from scipy.stats.kde import gaussian_kde
    CAN_PLOT_KDE = True
except ImportError:
    print "KDE will be omitted from window distribution plot due to missing scipy library."
    CAN_PLOT_KDE = False

def create_figure(figsize=(10,8)):
    fig = matplotlib.pyplot.figure(figsize=figsize)
    fig.hold = True
    return fig

def add_axis_to_figure(fig, subplot_layout=111, sharex=None, sharey=None):
    return fig.add_subplot(subplot_layout, sharex=sharex, sharey=sharey)

def save_figure(fig, fig_name):
    fig.savefig(fig_name)

def simple_plot(xs, ys, es, xlabel, ylabel):
    fig = create_figure()
    ax = add_axis_to_figure(fig)
    plot(ax, xs, ys, es, linewidth=1.0, xlabel=xlabel, ylabel=ylabel)

def plot(ax, xs, ys, es, xlabel=None, ylabel=None, label=None, linewidth=1, alpha=1, line_style="-", show=True):
    marker_size = 5
    width = 1.25 # table borders and ticks
    font_size = 12
    font_weight = "bold"

    line_color = "k"
    fill_color = "k"

    line_style = line_style

    if not ylabel:
        ylabel="y"
    if not xlabel:
        xlabel="x"
    if not label:
        label=""

    ax.errorbar(xs, ys, yerr=es, marker="o", zorder=1, linestyle=line_style, color=line_color, markerfacecolor=fill_color, label=label, markersize=marker_size, alpha=alpha, linewidth=linewidth)

    pylab.ylabel(ylabel, fontsize=font_size, fontweight=font_weight)
    pylab.xlabel(xlabel, fontsize=font_size, fontweight=font_weight)

    pylab.xticks(fontsize=font_size, fontweight=font_weight)
    pylab.yticks(fontsize=font_size, fontweight=font_weight)

    [i.set_linewidth(width) for i in ax.spines.itervalues()]

    if show:
        pylab.show()

def plot_with_kde(ax, xs, ys, xlabel=None, ylabel=None, label=None, show=True, kde=None, color="blue"):
    width = 1.25 # table borders and ticks
    font_size = 12
    font_weight = "bold"

    if not ylabel:
        ylabel="y"
    if not xlabel:
        xlabel="x"
    if not label:
        label=""

    ax.plot(xs, ys, color=color, label=label)
    if kde:
        kde_curve = kde(np.array(xs))
        ax.fill_between(xs, kde_curve, alpha=0.5, facecolor=color)

    pylab.ylabel(ylabel)
    pylab.xlabel(xlabel)

    pylab.xticks(fontsize=font_size, fontweight=font_weight)
    pylab.yticks(fontsize=font_size, fontweight=font_weight)

    [i.set_linewidth(width) for i in ax.spines.itervalues()]
    if show:
        pylab.show()

def plot_hist_to_axis(ax, centers, his, kde=None):
    plot_with_kde(ax, centers, his, show=False, kde=kde)

def finalise_figure(fig, ax, xlabel="x", ylabel="y", fig_name=None, show=False):
    pylab.ylabel(ylabel)
    pylab.xlabel(xlabel)

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

def plot_histogram(values, xlabel="x", ylabel="y", n_bins=100, fig_name=None, show=False):
    his, bins = np.histogram(values, bins = n_bins, normed=True)
    centers = (bins[:-1]+bins[1:])/2

    fig = create_figure((8,6))
    ax = add_axis_to_figure(fig)
    plot(ax, centers, his, np.zeros(len(his)), xlabel=xlabel, ylabel=ylabel, show=False)
    if fig_name:
        save_figure(fig, fig_name)
    if show:
        pylab.show()
    matplotlib.pyplot.close(fig)
