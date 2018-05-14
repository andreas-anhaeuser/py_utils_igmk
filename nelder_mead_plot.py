#!/usr/bin/python
"""Plot result of nelder_mead.minimize_nelder_mead

    Author
    ======
    Andreas Anhaeuser (AA)
    University of Cologne
    <andreas.anhaeuser@posteo.net>
"""

import numpy as np
from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt

_actionlabels = ('init', 'shr.', 'contr.', 'refl.', 'exp.')
_actionmarkers = ('d', '.', 'o', '^', 's')

def plot_all_parameters(record, filename, fun_log=False, fun_min=None,
        fun_max=None):
    best_params = record['best_vertices']
    best_fvals = record['best_fun_values']
    worst_params = record['worst_vertices']
    worst_fvals = record['worst_fun_values']
    actions = record['action']

    K, N = np.shape(best_params)
    Nplots = N + 1
    Nrows = int(np.ceil(np.sqrt(Nplots)))
    Ncols = int(np.ceil(1. * Nplots / Nrows))

    figsize = (297/25.4, 210/25.4)
    plt.figure(figsize=figsize)

    col_f = (0.5,) * 3
    col_param = 'r'
    col_action = (0.3, 0.3, 1.)
    lw_thick = 1
    lw_thin = 0.5
    lw_grid = 0.25
    ls_grid = ':'

    ms_action = 5

    fs_fun = 8

    if fun_min is None:
        if fun_log:
            fun_min = np.nanmin(best_fvals)
        else:
            fun_min = 0
    if fun_max is None:
        fun_max = np.nanmax(worst_fvals)

    for nplot in range(1, Nplots+1):
        ax1 = plt.subplot(Nrows, Ncols, nplot)
        
        ###################################################
        # FUNCTION VALUES                                 #
        ###################################################
        ax2 = ax1.twinx()
        if fun_log:
            # logarithmic
            plot_best = np.log10(best_fvals)
            plot_worst = np.log10(worst_fvals)
            ylabel = 'log10(function value)'
            ymin = np.log10(fun_min)
            ymax = np.log10(fun_max)
        else:
            # linear
            plot_best = best_fvals
            plot_worst = worst_fvals
            ylabel = 'function value'
            ymin = fun_min
            ymax = fun_max
        ax2.fill_between(np.arange(K), plot_best, ymin, color=col_f)
        ax2.plot(plot_worst, '--', color=col_f, lw=lw_thin)
        ax2.set_ylim(ymin, ymax)
        ax2.set_ylabel(ylabel, fontsize=fs_fun, color=col_f)
        ax2.tick_params('y', colors=col_f)
        # ax2.grid(color=col_f, lw=lw_grid/2, ls=ls_grid)
        ax2.set_zorder(0)

        if nplot == Nplots:
            color1 = col_action
            ax1.plot(actions, '-', lw=lw_thick, color=col_action)
            for n, action in enumerate(actions):
                ax1.plot(n, action, marker=_actionmarkers[action],
                        color=col_action, ms=ms_action)
            ax1.set_yticks(range(len(_actionlabels)))
            ax1.set_yticklabels(_actionlabels)
            plt.title('simplex transformation')
        else:
            color1 = col_param
            nparam = nplot - 1
            ax1.plot(best_params[:, nparam], '.-', color=col_param,
                    lw=lw_thick, zorder=11)
            ax1.plot(worst_params[:, nparam], '-', color=col_param,
                lw=lw_thin, zorder=10)
            plt.title('parameter %i' % nparam)

        ax1.tick_params('y', colors=color1)
        ax1.grid(color=color1, lw=lw_grid, ls=ls_grid)
        ax1.set_xlim(0, K-1)
        ax1.set_zorder(10)
        ax1.patch.set_visible(False)
        

    plt.tight_layout()
    plt.savefig(filename)
    print('Saved plot to %s' % filename)
    plt.close()

def plot_gradient_line(xx, yy, ff, vmin, vmax, cmap, N=3, lw=1):
    # interpolate
    x1, x2 = xx
    y1, y2 = yy
    f1, f2 = ff
    inc = 1./N
    a = np.arange(0, 1 + inc/2., inc)
    x = x1 + a * (x2 - x1)
    y = y1 + a * (y2 - y1)
    f = f1 + a * (f2 - f1)
    for n in range(len(a)-1):
        fmean = np.mean((f[n], f[n+1]))
        cval = (fmean - vmin) / (vmax - vmin)
        color = cmap(cval)
        plt.plot((x[n], x[n+1]), (y[n], y[n+1]), '-', color=color, lw=lw)


def plot_all_simplices(record, filename, axes=(1, 2)):
    a0, a1 = axes
    all_simplices = record['all_simplices']
    K, M, N = np.shape(all_simplices)
    fvals = record['all_fun_values']
    vmax = np.nanmax(record['worst_fun_values'])
    vmin = np.nanmin(record['best_fun_values'])
    vinc = (vmax - vmin) / 50.
    vlevels = np.arange(vmin, vmax + vinc / 2., vinc)
    cmap = plt.get_cmap('gist_rainbow', 2**10)

    figsize = (297/25.4, 210/25.4)
    plt.figure(figsize=figsize)

    for k in range(K):
        f = 0.9 * (1 - 1. * k/K)
        col_edge = 'k'

        x = all_simplices[k, :, a0]
        y = all_simplices[k, :, a1]
        f = fvals[k]

        for m1 in range(M):
            for m2 in range(m1 + 1, M):
                plot_gradient_line((x[m1], x[m2]), (y[m1], y[m2]), (f[m1], f[m2]),
                        vmin, vmax, cmap)
            plt.plot(x[m1], y[m1], 'k.')

    plt.xlabel('parameter %i' % a0)
    plt.ylabel('parameter %i' % a1)
    plt.grid()

    ###################################################
    # COLORBAR                                        #
    ###################################################
    ax = plt.gca()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    cb = plt.colorbar(sm)
    cb.set_label('function value')

    print('Saving to file %s ...' % filename)
    plt.savefig(filename)
    print('Done.')
    plt.close()


###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    filename = '/home/anhaeus/tmp/nmplot.eps'
    plot_all_parameters(record, filename, fun_log=True)
    filename = '/home/anhaeus/tmp/nm_simplices.eps'
    plot_all_simplices(record, filename)
