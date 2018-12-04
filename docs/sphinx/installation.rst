Installation
============

Installing Dependencies
-----------------------

Reaktoro has several software and library dependencies that need to be pre-installed for its successful compilation and installation. Dependencies varies across operating systems, which makes it a more painful task. However, we have integrated Reaktoro with a modern and powerful package dependency management system that simplifies the building process for Windows, Mac, and Linux operating systems. With this special addition to the project, all you have to install by yourself is `Conda`_!

Here are the steps for installing Conda:

1. Install `Miniconda <https://conda.io/miniconda.html>`_. We recommend a 64-bit installer for your operating system with Python version 3 (but if you prefer, an installer with Python version 2 should also work).
2. Add the conda channel **conda-forge**: ``conda config --append channels conda-forge``.
3. Install **conda-devenv**: ``conda install -n base conda-devenv``

Downloading Reaktoro
--------------------

We need now to download the source code of Reaktoro. This can be done by either issuing the following `git`_ command from the terminal:

.. code::

    git clone https://github.com/reaktoro/reaktoro.git

or by directly downloading the `latest version`_ in a zip file. If you use the direct download option, please unzip the downloaded file in a directory of your choice (assuming here the unzipped folder is named ``reaktoro`` for the next steps).

Creating and Activating a Conda Environment
-------------------------------------------

We need now to create and activate a specific conda environment that contains all the library dependencies needed to build Reaktoro.

1. In a terminal, change directory to the root of the source code: ``cd reaktoro``.
2. Create the conda environment: ``conda devenv``
3. Activate the environment: ``source activate reaktoro``.

.. note::

    If you are in Windows, you will actually need to execute ``activate reaktoro`` instead of ``source activate reaktoro``, which works only on Linux and MacOS.




Building and Installing Reaktoro
--------------------------------



.. code::

    cd reaktoro

You need now to setup



.. code::

    mkdir build && cd build
    cmake ..
    cmake --build .





Install Miniconda, pick the 64-bit installer that uses the latest Python version from: `conda.io/miniconda.html <https://conda.io/miniconda.html>`_.

.  The simplest way to install all required dependencies and successfully build Reaktoro is by using `Conda`.  from source get However, we have come up with a relatively simple way to build Reaktoro



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

Alternatively, you can access this `latest version`_ to directly download Reaktoro
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

.. _Conda: https://conda.io/docs/
.. _git: https://git-scm.com/
.. _latest version: https://github.com/reaktoro/reaktoro/archive/master.zip
.. _CVODE: https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode
.. _PHREEQC: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
.. _GEMS: http://gems.web.psi.ch/
.. _GitHub repository: https://github.com/reaktoro/reaktoro
.. _CMake: https://cmake.org/
.. _Gnuplot: http://www.gnuplot.info/