import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import os
from numpy import *

# Step 2: Auxiliary time related constants
second = 1
minute = 60
hour = 60 * minute
day = 24 * hour
year = 365 * day

# Step 10: Return a string for the title of a figure in the format Time: #h##m
def titlestr(t):
    t = t / minute   # Convert from seconds to minutes
    h = int(t) / 60  # The number of hours
    m = int(t) % 60  # The number of remaining minutes
    return 'Time: %2dh %2dm' % (h, m)

# Step 9: Generate figures for a result file
def plotfile_rt_kinetic_dolomitization(folder, file, x, nsteps, dt, output_quantities_map):

    step = int((file.split('.')[0]).split('-')[1])

    print('Plotting figure', step, '...')

    t = step * dt
    filearray = loadtxt(folder + '/' + file, skiprows=1)
    data = filearray.T

    ndigits = len(str(nsteps))

    plt.figure()
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=2.5, top=10.5)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('pH')
    plt.plot(x, data[output_quantities_map["pH"]])
    plt.tight_layout()
    plt.savefig('figures-' + folder + '/ph/{}.png'.format(str(step).zfill(ndigits)))
    plt.close()

    plt.figure()
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=40.0, top=50.0)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Mineral Volume [%$_{\mathsf{vol}}$]')
    plt.plot(x, data[output_quantities_map["phaseVolume(Calcite)"]] * 100, label='Calcite')
    plt.plot(x, data[output_quantities_map["phaseVolume(Dolomite)"]] * 100, label='Dolomite')
    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig('figures-' + folder + '/calcite-dolomite/{}.png'.format(str(step).zfill(ndigits)))
    plt.close()

    plt.figure()
    plt.yscale('log')
    plt.xlim(left=-0.02, right=0.52)
    plt.ylim(bottom=0.5e-5, top=2)
    plt.title(titlestr(t))
    plt.xlabel('Distance [m]')
    plt.ylabel('Concentration [molal]')
    plt.plot(x, data[output_quantities_map["speciesMolality(H+)"]], label='H+')
    plt.plot(x, data[output_quantities_map["speciesMolality(Ca++)"]], label='Ca++')
    plt.plot(x, data[output_quantities_map["speciesMolality(Mg++)"]], label='Mg++')
    plt.plot(x, data[output_quantities_map["speciesMolality(HCO3-)"]], label='HCO3-')
    plt.plot(x, data[output_quantities_map["speciesMolality(CO2(aq))"]], label='CO2(aq)')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('figures-' + folder + '/aqueous-species/{}.png'.format(str(step).zfill(ndigits)))
    plt.close()

    plt.close('all')

# Step 8: Plot all result files and generate a video
def plot_rt_kinetic_dolomitization(folder, x, nsteps, dt, output_quantities_map):

    ndigits = len(str(nsteps))

    # Plot all result files
    files = sorted(os.listdir(folder))

    Parallel(n_jobs=16)(delayed(plotfile_rt_kinetic_dolomitization)(folder, file, x, nsteps, dt, output_quantities_map) for file in files)
    # Create videos for the figures
    ffmpegstr = 'ffmpeg -y -r 30 -i figures-' + folder + '/{0}/%0' + str(ndigits) + 'd.png -codec:v mpeg4 -flags:v +qscale -global_quality:v 0 videos-' + folder + '/{0}.mp4'
    os.system(ffmpegstr.format('calcite-dolomite'))
    os.system(ffmpegstr.format('aqueous-species'))
    os.system(ffmpegstr.format('ph'))
