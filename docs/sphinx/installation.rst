Installation
============

Installing Dependencies
-----------------------

Reaktoro has several software and library dependencies that need to be
pre-installed for its successful compilation and installation. To greatly
simplify the building process of Reaktoro for Windows, macOS, and Linux
operating systems, we have integrated Reaktoro with a powerful package
dependency manager: `Conda`_.

Conda is the only dependency you should actually have to install manually to
build Reaktoro. The ``conda`` application is available by either installing
`Anaconda <https://www.anaconda.com/download>`_ or `Miniconda
<https://conda.io/miniconda.html>`_.

.. note::

    We recommend the installation of **Miniconda** instead of Anaconda, unless
    you think you need 1,400+ software packages that ship with Anaconda!
    Miniconda is just a tiny subset of Anaconda containing only our needed
    ``conda`` application and its dependencies. We also recommend the download
    of a **64-bit installer** corresponding to your operating system with
    **Python 3.7** or a newer version.

After installing Miniconda, as suggested above, go to a terminal and execute:

.. code:: bash

    conda config --append channels conda-forge
    conda install -n base conda-devenv

The first command above appends the `conda-forge <https://conda-forge.org/>`_
channel to the list of channels available in your ``conda`` installation. The
second command installs `conda-devenv <https://github.com/ESSS/conda-devenv>`_,
a conda tool with convenient functionalities to define and initialize conda
environments.


Downloading Reaktoro
--------------------

We need now to download the source code of Reaktoro, which is hosted at
`GitHub`_. This can be done by either executing the following `git`_ command
from the terminal (if you already have git installed!):

.. code:: bash

    git clone https://github.com/reaktoro/reaktoro.git

or by directly downloading the `latest version`_ of Reaktoro's source code in a
zip file.

.. note::

    If you use the direct download option above, please unzip the downloaded
    file in a directory of your choice. We assume the unzipped folder is named
    ``reaktoro`` for the next installation steps.


Creating and Activating the Reaktoro Conda Environment
------------------------------------------------------

The next step is to create and activate a conda environment that contains all
the software and library dependencies needed to build Reaktoro. In a terminal,
execute:

.. code:: bash

    cd reaktoro
    conda devenv
    source activate reaktoro

.. note::

    If you are in Windows, you will actually need to execute ``activate
    reaktoro`` instead of ``source activate reaktoro``, which works only in
    Linux and macOS.

The first command above changes directory to the root of the source code
directory. You might need to adapt the path to ``reaktoro``. The second command
initializes the conda environment with all software and library dependencies
needed to build Reaktoro. This can take a few minutes for the first time. The
third command activates the created *reaktoro environment*.

.. note::

    If you are curious about the list of dependencies needed to build Reaktoro,
    have a look at the file :download:`environment.devenv.yml
    <../../environment.devenv.yml>` in the root directory of Reaktoro's source
    code. This file is a *recipe* for the creation of a conda environment
    containing all dependencies needed by Reaktoro.


Building and Installing Reaktoro -- Easy Way
--------------------------------------------

Once the *Reaktoro conda environment* has been activated, as described above,
execute from the terminal:

.. code-block:: bash

    cmake -P install

.. note::

    This will first build and install Reaktoro in your *miniconda* directory if
    a conda environment is active, and not in a system directory such as
    ``/usr/local/`` in Linux. To install Reaktoro in a different place, say,
    under home directory, use: ``cmake -DPREFIX=$HOME -P install``.

.. attention::

    Compiling the Reaktoro C++ library and the Reaktoro Python module should
    take a few minutes.

Building and Installing Reaktoro -- Traditional Way
---------------------------------------------------

To build and install Reaktoro in a more traditional way, execute the following
from the root directory of Reaktoro's source code:

 .. code:: bash

    mkdir build && cd build
    cmake ..
    cmake --build . --target install


.. Install Miniconda, pick the 64-bit installer that uses the latest Python version from: `conda.io/miniconda.html <https://conda.io/miniconda.html>`_.

.. .  The simplest way to install all required dependencies and successfully build Reaktoro is by using `Conda`.  from source get However, we have come up with a relatively simple way to build Reaktoro



.. In the steps below we show how one can download Reaktoro, build, and
.. install it in Linux and macOS systems. Note that compiling Reaktoro can
.. take some time. This is because it heavily relies on template
.. metaprogramming for efficient vector and matrix calculations, as well as
.. for calculation of partial derivatives of most thermodynamic properties,
.. such as activity coefficients, phase molar volumes, standard Gibbs
.. energies, etc. In addition, it will also compile several third-party
.. libraries, such as `CVODE`_ for efficient solution of ordinary
.. differential equations (ODE), and the geochemical codes `PHREEQC`_ and
.. `GEMS`_. Compilation of the Python wrappers can also take several
.. minutes, as [Boost.Python] too relies on template metaprogramming.

.. Downloading the source code
.. ---------------------------

.. Reaktoro’s source code is kept in a `GitHub`_. If you have
.. ``git`` installed in your system, then downloading this repository is as
.. easy as running the following command in a terminal:

.. .. code:: bash

..    git clone https://github.com/reaktoro/reaktoro.git Reaktoro

.. Alternatively, you can access this `latest version`_ to directly download Reaktoro
.. source code as a zipped file. If you choose this option, unzip the file
.. before proceeding to the next step.

.. Installing the dependencies
.. ---------------------------

.. Reaktoro has a few dependencies that need to be installed before it can
.. be built. If you plan to compile only its C++ libraries, all you will
.. need is `CMake`_, which is used for managing and automating the whole
.. build process, including the installation of third party libraries. The
.. table below describes how to install CMake from the terminal in some
.. Linux distributions:

.. ========== ==============================
.. OS         Command
.. ========== ==============================
.. Ubuntu     ``sudo apt-get install cmake``
.. Fedora     ``sudo yum install cmake``
.. Arch Linux ``pacman -Ss cmake``
.. ========== ==============================

.. Optionally, you might want to install `Gnuplot`_ if you intend to do
.. real-time plotting of the chemical calculations.

.. ========== =============================================
.. OS         Command
.. ========== =============================================
.. Ubuntu     ``sudo apt-get install gnuplot5 gnuplot5-qt``
.. Fedora     ``sudo yum install gnuplot gnuplot-qt``
.. Arch Linux ``pacman -Ss gnuplot``
.. ========== =============================================

.. Check if a plot is successfuly output to a window terminal by issuing
.. the command:

.. .. code:: bash

..    gnuplot -persist -e 'plot sin(x)'

.. If a window did not show up with an interactive plot, you might need to
.. install a different package other than ``gnuplot-qt``. Check your
.. distribution, or install Gnuplot from source.

.. If you plan to use Reaktoro from Python, then a few more dependencies
.. are needed to compile the Python wrappers of Reaktoro’s C++ classes and
.. methods. If you just want the C++ libraries, you can skip this and

.. _Conda: https://conda.io/docs/
.. _git: https://git-scm.com/
.. _latest version: https://github.com/reaktoro/reaktoro/archive/master.zip
.. _CVODE: https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode
.. _PHREEQC: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
.. _GEMS: http://gems.web.psi.ch/
.. _GitHub: https://github.com/reaktoro/reaktoro
.. _CMake: https://cmake.org/
.. _Gnuplot: http://www.gnuplot.info/
