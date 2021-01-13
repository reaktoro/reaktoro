import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import font_manager as fm, rcParams
fpath = os.path.join(rcParams["datapath"], "texgyreadventor-regular.otf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(14)

import matplotlib as mpl
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 14
mpl.set_loglevel("critical")

# Plotting params
circ_area = 6 ** 2
custom_font = { }

C0 = '#107ab0'
C1 = '#fc5a50'

T_left = 25
T_right = 90
temperature_points = 14
temperatures = np.linspace(T_left, T_right, temperature_points)  # the x-coordinates of the plots

def empty_marker(color):
    return {'facecolor': 'white', 'edgecolor': color, 's': circ_area, 'zorder': 2, 'linewidths': 1.5 }

def filled_marker(color):
    return {'color': color, 'marker': 'D', 'zorder': 2}

def line_empty_marker(color):
    return {'marker': 'd', 'markerfacecolor': 'white', 'markeredgecolor':color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line_filled_marker(color):
    return {'color': color, 'markersize': 6, 'markeredgewidth': 1.5 }

def line(color):
    return {'linestyle': '-', 'color': color, 'zorder': 1, 'linewidth': 2}

def line_error(color):
    return {'linestyle': ':', 'marker': 'D', 'color': color, 'zorder': 1, 'linewidth': 0.5, 'markersize': 4}

def plot_concentrations(data_reaktoro, data_phreeqc, molalities, tag, title, folder):

    colors = ['C1', 'C2', 'C3']
    plt.axes(xlim=(temperatures[0] - 2, temperatures[-1] + 2))

    for molality, i in zip(molalities, list(range(len(temperatures)))):

        plt.xlabel(r'Temperature [°C]', fontproperties=prop)
        plt.ylabel('Solubility [mol/kgw]', fontproperties=prop)
        plt.title(title, fontproperties=prop)

        plt.plot(temperatures, data_reaktoro[i], label='Reaktoro, ' + str(molality) + ' molal', **line(colors[i]))
        plt.plot(temperatures, data_phreeqc[i], 'D', **line_filled_marker(colors[i]))[0],
        plt.plot([], [], 'D', label='PHREEQC, ' + str(molality) + ' molal', **line_filled_marker(colors[i]))

    plt.legend(loc='upper right', prop=prop)
    plt.savefig(folder + '/reaktoro-phreeqc-comparison' + tag +'.pdf')
    plt.close()

if __name__ == '__main__':

    molalities = [1.0, 2.0, 4.0]

    results_folder = 'results-kineticpath-complex-souring'
    results_file = results_folder + '/kineticpath-scavenging-complex-t0-0-tfinal-5040000-n-1400-dh-full-conv-kin-conv-eq.txt'

    # Load the Reaktoro values from the file
    data_reaktoro_pyrrhotite = np.loadtxt(results_file, skiprows=2, usecols=(-6))
    data_reaktoro_calcite = np.loadtxt(results_file, skiprows=2, usecols=(-5))
    data_reaktoro_quartz = np.loadtxt(results_file, skiprows=2, usecols=(-4))
    data_reaktoro_kaolinite = np.loadtxt(results_file, skiprows=2, usecols=(-3))
    data_reaktoro_daphnite = np.loadtxt(results_file, skiprows=2, usecols=(-2))
    data_reaktoro_siderite = np.loadtxt(results_file, skiprows=2, usecols=(-1))
    data_reaktoro_time = np.loadtxt(results_file, skiprows=2, usecols=(0))

    # Load the PHREEQC values from the file -d_CO2(g) concentrations
    results_file = results_folder + '/H2S_scavenging_with_kinetic_calcite_1400.txt'
    data_phreeqc_pyrrhotite = np.loadtxt(results_file, skiprows=2, usecols=(-22))
    data_phreeqc_calcite = np.loadtxt(results_file, skiprows=2, usecols=(-11))
    data_phreeqc_quartz = np.loadtxt(results_file, skiprows=2, usecols=(-9))
    data_phreeqc_kaolinite = np.loadtxt(results_file, skiprows=2, usecols=(-7))
    data_phreeqc_daphnite = np.loadtxt(results_file, skiprows=2, usecols=(-5))
    data_phreeqc_siderite = np.loadtxt(results_file, skiprows=2, usecols=(-3))

    #print(data_phreeqc)
    #input()

    tag = ''
    nsteps = len(data_reaktoro_time)
    step = 20

    # ---------------------------------------------------------------------------------------------------------------- #
    # Plot siderite and pyrrhotite
    # ---------------------------------------------------------------------------------------------------------------- #
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_xlim([data_reaktoro_time[0] - 1e5, data_reaktoro_time[-1] + 1e5])
    ax1.set_ylim([0.0875, 0.101])
    ax1.set_ylabel('Amount [mol]', fontproperties=prop)
    ax1.set_xlabel('Time [s]', fontproperties=prop)
    ax1.plot(data_reaktoro_time[0:nsteps:step], data_reaktoro_siderite[0:nsteps:step], label='Siderite (Reaktoro)', color="C2")
    ax1.plot(data_reaktoro_time[0:nsteps:step], data_phreeqc_siderite[0:nsteps:step], 'o', label='Siderite (PHREEQC)', **filled_marker('C2'))
    ax1.legend(loc='best', prop=prop)

    ax2.set_xlim([data_reaktoro_time[0] - 1e5, data_reaktoro_time[-1] + 1e5])
    ax2.set_ylim([-0.001, 0.0125])
    ax2.set_ylabel('Amount [mol]', fontproperties=prop)
    ax2.plot(data_reaktoro_time[0:nsteps:step], data_reaktoro_pyrrhotite[0:nsteps:step], label='Pyrrhotite (Reaktoro)', color="C1")
    ax2.plot(data_reaktoro_time[0:nsteps:step], data_phreeqc_pyrrhotite[0:nsteps:step], 'o', label='Pyrrhotite (PHREEQC)', **filled_marker('C1'))
    ax2.legend(loc='best', prop=prop)

    plt.xlabel('Time [s]', fontproperties=prop)
    plt.savefig(results_folder + '/reaktoro-phreeqc-pyrrhotite-siderite-comparison.png')
    #plt.show()
    plt.close()

    # ---------------------------------------------------------------------------------------------------------------- #
    # Plot daphnite and kaolinite
    # ---------------------------------------------------------------------------------------------------------------- #
    fig = plt.figure()

    plt.axes(xlim=(data_reaktoro_time[0] - 1e5, data_reaktoro_time[-1] + 1e5))
    plt.xlabel('Time', fontproperties=prop)
    plt.ylabel('Amount [mol]', fontproperties=prop)
    plt.plot(data_reaktoro_time[0:nsteps:step], data_reaktoro_daphnite[0:nsteps:step], label='Daphnite (Reaktoro)', color="C3")
    plt.plot(data_reaktoro_time[0:nsteps:step], data_reaktoro_kaolinite[0:nsteps:step], label='Kaolinite (Reaktoro)', color="C4")
    #plt.plot(data_reaktoro_time[0:nsteps:step], data_reaktoro_calcite[0:nsteps:step], label='Calcite (Reaktoro)', color="C6")

    plt.plot(data_reaktoro_time[0:nsteps:step], data_phreeqc_daphnite[0:nsteps:step], 'o', label='Daphnite (PHREEQC)', **filled_marker('C3'))
    plt.plot(data_reaktoro_time[0:nsteps:step], data_phreeqc_kaolinite[0:nsteps:step], 'o', label='Kaolinite (PHREEQC)', **filled_marker('C4'))
    #plt.plot(data_reaktoro_time[0:nsteps:step], data_phreeqc_calcite[0:nsteps:step], 'o', label='Calcite (PHREEQC)', **filled_marker('C6'))

    plt.legend(loc='best', prop=prop)
    #plt.show()
    plt.savefig(results_folder + '/reaktoro-phreeqc-daphnite-kaolinite-comparison.png')
    plt.close()

    # ---------------------------------------------------------------------------------------------------------------- #
    # Plot calcite
    # ---------------------------------------------------------------------------------------------------------------- #
    fig = plt.figure()

    plt.axes(xlim=(data_reaktoro_time[0] - 1e5, data_reaktoro_time[-1] + 1e5))
    plt.xlabel('Time', fontproperties=prop)
    plt.ylabel('Amount [mol]', fontproperties=prop)
    plt.plot(data_reaktoro_time[0:nsteps:step], data_reaktoro_calcite[0:nsteps:step], label='Calcite (Reaktoro)', color="C6")
    plt.plot(data_reaktoro_time[0:nsteps:step], data_phreeqc_calcite[0:nsteps:step], 'o', label='Calcite (PHREEQC)', **filled_marker('C6'))

    plt.legend(loc='best', prop=prop)
    #plt.show()
    plt.savefig(results_folder + '/reaktoro-phreeqc-calcite-comparison.png')
    plt.close()