Compiling on Microsoft Visual Studio Community 2013
===================================================

Necessary Changes to PHREEQC
----------------------------

Unless `struct phase` in file `global_structures.h` is not changed to a class, `class phase`, with public interface, MSVC resolves to a LNK2019 unresolved symbol error on method:

~~~
LDBLE calc_PR(std::vector<class phase *> phase_ptrs, LDBLE P, LDBLE TK, LDBLE V_m);
~~~

where `class phase` should in fact be `struct phase`. Thus, changing the definition of `phase` fixes this issue.

Python and Numpy Installation
-----------------------------

Download the most recent version of [Python](https://www.python.org/downloads/windows/). If the x86-64 was downloaded, then a x64 version of Numpy must be installed, which can be found [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/). Choose the one compiled with MKL.

After the installation of Python, which also probably installed `pip`, go to the directory where numpy was downloaded. Then issue the command:

~~~
pip install numpy‑1.9.2+mkl‑cp27‑none‑win_amd64.whl
~~~

where the example above is for the numpy binaries for Python v2.7, compiled with MKL, for Windows x64.

Boost Installation
------------------

Boost in Windows can be installed using the binaries provided [here](http://sourceforge.net/projects/boost/files/boost-binaries/). Ensure you choose the most recent version and the appropriate file for your OS (x32 or x64).

CMake Errors About Boost
------------------------

When setting up the build configuration of Reaktoro for MSVC 2013, the Boost libraries might not be found automatically by CMake. In this case, provide the path to the Boost directory using the definitions BOOST_ROOT. You will also need to provide the path to the directory where the compiled boost libraries are (e.g., libboost_python). This will require setting the BOOST_LIBRARYDIR. For example, is Boost is installed under `C:\local\boost_1_58_0` and the libraries are found in `C:\local\boost_1_58_0\lib64-msvc-12.0`, then set these paths to BOOST_ROOT and BOOST_LIBRARYDIR respectively.

If you get some python error during CMake configure step, ensure (1) the Python directory was set in the environmental variable `Path` of Windows, (2) numpy has been installed, (3) the installed Python has the same architecture as Boost (x64 or x32).


Download Boost binaries compatible with the MSVC version installed in the system.

Install gnuwin32, nsis.

C:\Python27-w32\Scripts\pip.exe install tabulate pyyaml
C:\Python27-w64\Scripts\pip.exe install tabulate pyyaml

set PATH=C:\gnuwin32\bin;C:\Python27-w32\;C:\Python27-w32\Scripts\;C:\boost32\lib32-msvc-12.0;%PATH%

cmake ../.. -DBUILD_ALL=ON -DBOOST_INCLUDEDIR=C:\boost32 -DBOOST_LIBRARYDIR=C:\boost32\lib32-msvc-12.0
cmake ../.. -DBUILD_ALL=ON -DBOOST_ROOT=C:\boost32

cmake --build . --config Release
cmake --build . --config Release --target package
