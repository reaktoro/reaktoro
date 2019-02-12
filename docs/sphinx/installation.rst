Installation
============

There are many ways for installing Reaktoro, from a simple single-line terminal
command to building the entire project from sources. All current possibilities
are detailed below.

Installation using Conda
------------------------

Reaktoro can be easily installed using `Conda`_, a powerfull package manager
used to simplify Reaktoro's installation and the management of its
external software dependencies. Once you install Conda, and append the
necessary channels, you'll be able to install Reaktoro by just executing the
following command in your terminal:

.. code::

    conda install reaktoro


.. _Anaconda: https://www.anaconda.com/download
.. _Miniconda: https://conda.io/miniconda.html


Follow the Conda installation steps shown next before you execute this command!

Installing Conda
^^^^^^^^^^^^^^^^

Conda can be installed by installing either Anaconda_ or Miniconda_. We
recommend the installation of Miniconda, unless you already have Anaconda
installed or you think you need 1,400+ software packages that ship with it!
Miniconda is just a tiny subset of Anaconda containing only our needed
``conda`` application and its dependencies.

Download the **Miniconda Installer** by clicking on the image below:

.. figure:: img/logos/conda-logo.svg
    :figwidth: 80%
    :width: 80%
    :align: center
    :target: Miniconda_

    Reaktoro is a proud user of Conda, a powerfull and modern package manager!


You'll be given options to install `Miniconda`_ for Windows, Mac OS X, and
Linux (32-bit or 64-bit) using either Python 3.7 or Python 2.7. We recommend a
**64-bit installer** of Miniconda with **Python 3.7**.

Adding conda-forge channels
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Reaktoro pre-built package is hosted on `conda-forge
<https://anaconda.org/conda-forge/reaktoro>`_. After installing Miniconda, go
to a terminal and execute:

.. code::

    conda config --append channels conda-forge
    conda config --append channels conda-forge/label/gcc7

to add the conda-forge channels required to find the Reaktoro package.

All should now be set to install Reaktoro using:

.. code::

    conda install reaktoro

Go to `Reporting a failed installation`_ if this does not work for you.

Installation using CMake
------------------------

Reaktoro has several software and library dependencies that need to be
pre-installed for its successful compilation and installation using CMake_. To
greatly simplify the building process of Reaktoro for Windows, Mac OS X, and
Linux, you'll need `Conda`_. Follow the Conda installation steps in the
previous section, in which a Miniconda installer is used.

After installing Miniconda, go to a terminal and execute:

.. code::

    conda install -n base conda-devenv

This installs `conda-devenv`_, a conda tool with convenient functionalities to
define and initialize conda environments.

Downloading Reaktoro from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We need now to download the source code of Reaktoro, which is hosted on
`GitHub`_. This can be done by either executing the following `git`_ command
from the terminal (if you already have git installed!):

.. code:: bash

    git clone https://github.com/reaktoro/reaktoro.git

.. |master.zip| replace:: :download:`reaktoro-master.zip<https://github.com/reaktoro/reaktoro/archive/master.zip>`

or by directly downloading |master.zip|, the latest version of Reaktoro's
source code in a zip file.


.. note::

    If you use the direct download option above, please unzip the downloaded
    file in a directory of your choice. We assume the unzipped folder is named
    ``reaktoro`` for the next installation steps, and not ``reaktoro-master``.

Creating a conda environment for Reaktoro
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next step is to create a **conda environment** that contains all
the software and library dependencies needed to build Reaktoro. In the root of
the reaktoro directory, execute:

.. code:: bash

    conda devenv

This command will create the conda environment called ``reaktoro``, which can
take a few minutes to complete for the first time.

.. attention::

    You only need to execute ``conda devenv`` again when the list of external
    dependencies changes or some configuration in the conda environment
    ``reaktoro`` is altered.

.. note::

    If you are curious about the list of dependencies needed to build Reaktoro,
    have a look at the file :download:`environment.devenv.yml
    <../../environment.devenv.yml>` in the root directory of Reaktoro's source
    code. This file is a *recipe* for the creation of our conda environment
    ``reaktoro`` containing all required dependencies.

Activating the conda environment for Reaktoro
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next step is to **activate the conda environment** ``reaktoro`` that
conda-devenv_ created for us:

.. code::

    conda activate reaktoro

.. attention::

    You need to activate the ``reaktoro`` conda environment whenever you use
    Reaktoro from C++ or Python! This is because conda will adjust some
    environment variables in your system (e.g., ``PYTHONPATH``,
    ``LD_LIBRARY_PATH``, ``PATH``) so that Reaktoro's libraries, executables,
    and Python packages can be found. Activating the ``reaktoro`` conda
    environment is the simplest way to get these environment variables set
    correctly.

Building and installing Reaktoro with CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can now build and install Reaktoro by executing the following from the root
of the reaktoro source directory:

.. code::

    cmake -P install

Assuming the conda environment ``reaktoro`` is active, this command will first
build Reaktoro and then install its header files, libraries, executables,
Python package in your local *miniconda* directory:

:Linux: ``/home/user/miniconda3/envs/reaktoro/``
:Windows: ``C:\miniconda3\envs\reaktoro\``

Alternatively, to build and install Reaktoro in a more traditional way, execute
the following from the root directory of Reaktoro's source code:

.. code::

    mkdir build
    cd build
    cmake ..
    cmake --build . --target install

The following is also possible with CMake v3.13 or newer:

.. code::

    cmake -S . -B build
    cmake --build build/ --target install

.. tip::

    Compiling the Reaktoro C++ library and the Reaktoro Python module should
    take a few minutes for the first time. However, if you activate the
    reaktoro conda environment, `ccache <https://ccache.samba.org/>`_ will be
    used to significantly speed up future compilations automatically for you!

Installing Reaktoro in a custom directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install Reaktoro in a different directory, say, ``/home/user/other``, use:

.. code::

    cmake -DPREFIX=/home/user/other -P install

or

.. code:: bash

    cmake .. -DCMAKE_INSTALL_PREFIX=/home/user/other
    cmake --build . --target install

You'll need, however, to set the environment variables ``PYTHONPATH``,
``LD_LIBRARY_PATH``, and ``PATH`` yourself. For example, in Linux:

.. code-block:: bash

    export PATH=$PATH:/home/user/other/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/other/lib{64}
    export PYTHONPATH=$PYTHONPATH:/home/user/other/lib{64}/pythonX.Y/site-packages

where ``lib{64}`` is either ``lib`` or ``lib64``, and ``pythonX.Y`` where
``pythonX.Y`` depends on the python version used to compile Reaktoro's Python
package (e.g., ``python3.6``, ``python3.7``).

Checking for a successful installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check if Reaktoro was installed correctly by executing:

.. code:: bash

    python -c 'import reaktoro; print(reaktoro.__path__[0])'

This should print the path to the installed python package ``reaktoro``. For
example:

.. code:: bash

    /home/user/miniconda3/envs/reaktoro/lib/pythonX.Y/site-packages/reaktoro

where ``pythonX.Y`` above depends on the python version used.

.. attention::

    Make sure you have the conda environment ``reaktoro`` active! Otherwise the
    checking above might not work without further actions (e.g., changing the
    ``PYTHONPATH`` environment variable).

If you get instead something like:

.. code:: bash

    Traceback (most recent call last):
        File "<string>", line 1, in <module>
    ModuleNotFoundError: No module named 'reaktoro'

then the installation was not successful or it was installed in a custom path
that is not yet given in the ``PYTHONPATH`` environment variable.

Reporting a failed installation
-------------------------------

A failed installation can be a frustating emotion, but we will be happy to help
you fixing your installation issue. However, please **do make sure** you
followed exactly the steps given before. If you are sure that you followed
every single instruction and the installation still fails, please go to:

.. centered::
    `Reaktoro's GitHub Issues`_

and let us know!

.. _CMake: https://cmake.org/
.. _conda-devenv: https://github.com/ESSS/conda-devenv
.. _Conda: https://conda.io/docs/
.. _git: https://git-scm.com/
.. _GitHub: https://github.com/reaktoro/reaktoro
.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/Reaktoro/issues/new
