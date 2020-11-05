"""Module for plotting variables against time or position"""

import matplotlib.pyplot as plt
import os
import seaborn as sns
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import pycott_not_link.cott_utils

path = os.path.abspath(__file__) # Figures out the absolute path, in case working directory moves

def plot_var_vs_t(varstr, icell, tec):
    """ Plot variable in given cell against time.

    Arguments
        varstr      Name of variable in tec-object.
        icell       Index of cell to get data from.
        tec         Knuplot Tecplot object
    """

    t, var = pycott.cott_utils.get_var_vs_t(varstr, icell, tec)

    plt.plot(t, var, ".-")
    plt.xlabel(r"$t$ (s)")
    plt.ylabel(format_text(varstr))
    plt.grid()
    plt.show()

    
def plot_vars_vs_x(varstr_list, t_list, tec, filename="var_vs_x.pdf",
                   show_plot=True, save_plot=False):
    """ Plot variables at given times against x.
    
    Arguments
        varstr_list    Name of variable(s) in tec-object (list of str)
        t_list         The time(s) variables should be plotted over x (list of doubles)
        tec       Knuplot Tecplot object
    """
   
    folder = 'plots/'
    filename = folder + filename
    
    x, var_vals = pycott.cott_utils.get_vars_vs_x(varstr_list, t_list, tec)

    n_vars = len(varstr_list)
    n_times = len(t_list)

    # Markers for plotting, can only plot 8 different times
    # in same plot. More markers can be added if needed.
    markers = ['-', '--', '-.', ':', '.', 'o', '+', 'x']

    fig, axs = plt.subplots(1, n_vars)
    if n_vars == 1:
        axs = [axs]
    for i in range(n_vars):
        for j in range(n_times):
            axs[i].plot(x, var_vals[j, i, :], markers[j],
                        label = ("t = " + str(t_list[j])+ " s"))
        axs[i].set_xlabel('x [m]')
        axs[i].set_ylabel(format_text(varstr_list[i]))
        axs[i].legend()

    
    plt.subplots_adjust(wspace=0.35)
    plt.tight_layout
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(filename)
    for ax in axs:
        ax.clear()
    plt.close()

    
def compare_var_vs_x(varstr, t_list, tec_list, names='', filename='comparison_over_x.pdf',
                     show_plot=True, save_plot=False):
    """ Plot one variable at given times against x for different simulations.
    Shows and/or saves a figure with a subplot for each time in the t_list,
    where the variable is plotted over x for the different simulations recorded in tec_list.
    
    Arguments
        varstr_list    Name of variable(s) in tec-object (list of str)
        t_list         The time(s) variables should be plotted over x (list of doubles)
        tec_list       Knuplot Tecplot object list
    """
   
    folder = 'plots/'
    filename = folder + filename

    n_sims = len(tec_list)
    n_times = len(t_list)

    # Markers for plotting, can only plot 8 different simulations in same plot
    markers = ['-', '--', '-.', ':', '.', 'o', '+', 'x']

    fig, axs = plt.subplots(1, n_times)
    if n_times == 1:
        axs = [axs]

    if names == '':
        names = []
        for sim in range(n_sims):
            names.append('simulation ' + str(sim))
            
    for i in range(n_times):
        for j in range(n_sims):
            x, var = pycott.cott_utils.get_vars_vs_x([varstr], [t_list[i]], tec_list[j])
            axs[i].plot(x, var[0, 0, :], markers[j],
                        label = names[j])
            
        axs[i].set_xlabel('x [m]')
        axs[i].set_ylabel(format_text(varstr))
        axs[i].set_title(("t = " + str(t_list[i])+ " s"))
        axs[i].legend()

    
    plt.subplots_adjust(wspace=0.35)
    plt.tight_layout
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(filename)
    for ax in axs:
        ax.clear()
    plt.close()

def compare_var_vs_t(varstr, pos_list, tec_list, names='', filename="comparison_over_t",
                     show_plot=True, save_plot=False):
    """ Plot one variable at given times against x for different simulations.
    Shows and/or saves a figure with a subplot for each time in the t_list,
    where the variable is plotted over x for the different simulations recorded in tec_list.
    
    Arguments
        varstr_list    Name of variable(s) in tec-object (list of str)
        pos_list         The position(s)[m] variables should be plotted over t (list of int?)
        tec_list       Knuplot Tecplot object list
        names          Possible list over names of the simulations eg ['TFM', 'HEM']
    """
   
    folder = 'plots/'
    filename = folder + filename

    n_sims = len(tec_list)
    n_pos = len(pos_list)

    # Markers for plotting, can only plot 8 different simulations in same plot
    markers = ['-', '--', '-.', ':', '.', 'o', '+', 'x']

    if names == '':
        names = []
        for sim in range(n_sims):
            names.append('simulation ' + str(sim))
            
    fig, axs = plt.subplots(1, n_pos)
    if n_pos == 1:
        axs = [axs]


    for i in range(n_pos):
        for j in range(n_sims):
            t, var, true_x = pycott.cott_utils.get_var_vs_t_at_x(varstr, pos_list[i],
                                                                 tec_list[j])
            axs[i].plot(t, var, markers[j], label = (names[j] + ' at x = ' + str(true_x)))
            
        axs[i].set_xlabel('t [s]')
        axs[i].set_ylabel(format_text(varstr))
        axs[i].set_title((r"x $\approx$ " + str(pos_list[i])))
        axs[i].legend()

    
    plt.subplots_adjust(wspace=0.35)
    plt.tight_layout
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(filename)
    for ax in axs:
        ax.clear()
    plt.close()

    
def plot_var_3D(varstr, tec, filename="var_xt_3Dplot", save_plot=False, show_plot=True):
    """
    Plots a variable over time and space in a 3D plot. Giving a surface plot.
    """
    folder = 'plots/'
    filename = folder + filename
    
    t, var = pycott.cott_utils.get_var_vs_t(varstr, 0, tec)
    x, var = pycott.cott_utils.get_vars_vs_x([varstr], [0], tec)
    # only do this to get x and t to create a meshgrid

    nx = len(x)

    x_grid, t_grid = np.meshgrid(x, t, sparse=False, indexing='xy')
    var3D = np.zeros_like(x_grid)
    for i in range(nx):
        temp, var3D[:, i] = pycott.cott_utils.get_var_vs_t(varstr, i, tec)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(x_grid, t_grid, var3D, cmap=plt.cm.viridis, linewidth=0.2)
    ax.set_xlabel(r'$x [m]$')
    ax.set_ylabel(r'$t [s]$')
    ax.set_title(format_text(varstr))

    plt.tight_layout()
    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig(filename)
    ax.clear()
    plt.close()


def format_text(varstr):
    if len(varstr) == 1:
        return ('$' + varstr + '$')
    else:
        formatted = ''
        for i in range(len(varstr)):
            if varstr[i] == '`':
                if varstr[i+1] ==  'r':
                    formatted += r'$\rho$'
                elif varstr[i+1] == 'a':
                    formatted += r'$\alpha$'
                elif varstr[i+1] == 'm':
                    formatted += r'$\mu$'
                elif formatted[i+1] == 's':
                    formatted += r'$\sigma$'
                elif varstr[i+1] == 'S':
                    formatted += r'$\Sigma$'
            elif varstr[i] == '_':
                if varstr[i+1] == str(1):
                    formatted += r'$_1$'
                elif varstr[i+1] == str(2):
                    formatted += r'$_2$'
                elif varstr[i+1] == str(3):
                    formatted += r'$_3$'
                elif varstr[i+1] == 'g':
                    formatted += r'$_g$'
                elif varstr[i+1] == 't':
                    formatted += r'$_t$'
                elif varstr[i+1] == 'i':
                    formatted += r'$_i$'
            elif varstr[i] == '^' and varstr[i+1] == '~':
                formatted = formatted[0: -1]
                if varstr[i-1] == 'c':
                    formatted += r'$\tilde{c}$'
            elif varstr[i-1] == '`':
                pass
            elif varstr[i-1] == '_':
                pass
            elif varstr[i-1] == '^':
                pass
            else:
                formatted += varstr[i]
        return formatted

    
if __name__ == "__main__":
    # Example:
    import sys
    varname = sys.argv[1]
    ic = int(sys.argv[2])
    tecpath = sys.argv[3]
    tecobj = pycott.cott_utils.load_tec(tecpath)
    plot_var_vs_t(varname, ic, tecobj)
    
