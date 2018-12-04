Installation
============

In the steps below we show how one can download Reaktoro, build, and
install it in Linux and MacOS systems. Note that compiling Reaktoro can
take some time. This is because it heavily relies on template
metaprogramming for efficient vector and matrix calculations, as well as
for calculation of partial derivatives of most thermodynamic properties,
such as activity coefficients, phase molar volumes, standard Gibbs
energies, etc. In addition, it will also compile several third-party
libraries, such as `CVODE`_ for efficient solution of ordinary
differential equations (ODE), and the geochemical codes `PHREEQC`_ and
`GEMS`_. Compilation of the Python wrappers can also take several
minutes, as [Boost.Python] too relies on template metaprogramming.

Downloading the source code
---------------------------

Reaktoro’s source code is kept in a `GitHub repository`_. If you have
``git`` installed in your system, then downloading this repository is as
easy as running the following command in a terminal:

.. code::

   git clone https://github.com/reaktoro/reaktoro.git Reaktoro

Alternatively, you can access this `link`_ to directly download Reaktoro
source code as a zipped file. If you choose this option, unzip the file
before proceeding to the next step.

Installing the dependencies
---------------------------

Reaktoro has a few dependencies that need to be installed before it can
be built. If you plan to compile only its C++ libraries, all you will
need is `CMake`_, which is used for managing and automating the whole
build process, including the installation of third party libraries. The
table below describes how to install CMake from the terminal in some
Linux distributions:

========== ==============================
OS         Command
========== ==============================
Ubuntu     ``sudo apt-get install cmake``
Fedora     ``sudo yum install cmake``
Arch Linux ``pacman -Ss cmake``
========== ==============================

Optionally, you might want to install `Gnuplot`_ if you intend to do
real-time plotting of the chemical calculations.

========== =============================================
OS         Command
========== =============================================
Ubuntu     ``sudo apt-get install gnuplot5 gnuplot5-qt``
Fedora     ``sudo yum install gnuplot gnuplot-qt``
Arch Linux ``pacman -Ss gnuplot``
========== =============================================

Check if a plot is successfuly output to a window terminal by issuing
the command:

.. code::

   gnuplot -persist -e 'plot sin(x)'

If a window did not show up with an interactive plot, you might need to
install a different package other than ``gnuplot-qt``. Check your
distribution, or install Gnuplot from source.

If you plan to use Reaktoro from Python, then a few more dependencies
are needed to compile the Python wrappers of Reaktoro’s C++ classes and
methods. If you just want the C++ libraries, you can skip this and

.. _CVODE: https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode
.. _PHREEQC: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
.. _GEMS: http://gems.web.psi.ch/
.. _GitHub repository: https://github.com/reaktoro/reaktoro
.. _link: https://github.com/reaktoro/reaktoro/get/master.zip
.. _CMake: https://cmake.org/
.. _Gnuplot: http://www.gnuplot.info/