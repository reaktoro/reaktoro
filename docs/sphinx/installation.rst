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

    If your conda version is not 4.4 or newer, and you are using Windows, then
    you will actually need to execute ``activate reaktoro`` instead of ``source
    activate reaktoro``.

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

The easy way shown before is more or less these sequence of commands!

Checking for a Successful Installation
--------------------------------------

Check if Reaktoro was installed correctly by executing:

.. code:: bash

    python -c 'import reaktoro; print(reaktoro.__file__)'.

This should print the path to the installed python package ``reaktoro``. For
example:

.. code:: bash

    /home/your-user-name/miniconda3/envs/reaktoro/lib/pythonX.Y/site-packages/reaktoro/__init__.py

where ``pythonX.Y`` above depends on the python version used.

.. attention::

    Make sure you have the **reaktoro conda environment active**! Otherwise the
    checking above might not work without further actions (e.g., changing
    the ``PYTHONPATH`` environment variable).

If you get instead something like:

.. code:: bash

    Traceback (most recent call last):
        File "<string>", line 1, in <module>
    ModuleNotFoundError: No module named 'reaktoro'

then the installation was not successful or it was installed in a custom path
that is not yet given in the ``PYTHONPATH`` environment variable.

My Installation Failed. What do I do?
-------------------------------------

This can be frustating, but we will be happy to help you fixing your
installation issue. However, please **do make sure** you followed exactly the
steps given previously. If you are sure that you followed every single
instruction and the installation still fails, please go to `Reaktoro's GitHub
Issues`_ and let us know (please give some details and describe your operating
system).


.. _Conda: https://conda.io/docs/
.. _git: https://git-scm.com/
.. _latest version: https://github.com/reaktoro/reaktoro/archive/master.zip
.. _GitHub: https://github.com/reaktoro/reaktoro
.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/Reaktoro/issues/new
