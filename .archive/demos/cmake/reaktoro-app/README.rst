Using CMake to Build a C++ Reaktoro App
=======================================

This example demonstrates how your C++ Reaktoro application (executable) can be
built using `CMake <https://cmake.org/>`_.

In the root directory of this example, execute the following in the terminal:
::

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

.. attention::

    If you have a custom installation of Reaktoro (i.e., Reaktoro is not in a
    path that CMake can find it), you'll need to use the following instead of
    just ``cmake ..`` above:
    ::

        cmake .. -DCMAKE_PREFIX_PATH=/home/username/reaktoro/

    assuming Reaktoro was installed in your home directory (``username`` needs
    to be changed!)

After building the application, you can run the app:

::

    $ ./app

which should output a file named ``state-cmake-example.txt``.

.. attention::

    If you are using Conda to manage the dependencies of Reaktoro, then you
    should activate the Reaktoro conda environment:
    ::

        conda activate reaktoro

    Read more on the installation instructions of Reaktoro.
