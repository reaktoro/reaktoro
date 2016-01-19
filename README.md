#Reaktoro

Reaktoro is a unified framework for modeling chemically reactive systems. It provides methods for chemical equilibrium and kinetics calculations for multiphase systems. Reaktoro is mainly developed in C++ for performance reasons. A Python interface is available for a more convenient and simpler use of the scientific library. Currently Reaktoro can interface with two widely used geochemical software: [PHREEQC](http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/) and [GEMS](http://gems.web.psi.ch/). This document describes how to download and install Reaktoro, and demonstrate some basic usage.

----------

##Installation

In the steps below we will show how one can download Reaktoro, build, and install it in Linux and MacOS systems. We plan to release binaries (i.e., the libraries already compiled) for Windows soon. Please get in touch so we can know how urgent these binaries are for you.

### Downloading Reaktoro
Reaktoro source code is kept in this [Bitbucket repository](https://bitbucket.org/reaktoro/reaktoro). If you have `git` installed in your system, then downloading this repository is as easy as running the following command in a terminal:

    git clone https://bitbucket.org/reaktoro/reaktoro.git Reaktoro

Alternatively, you can access this [link](https://bitbucket.org/reaktoro/reaktoro/get/master.zip) to directly download Reaktoro source code as a zipped file. If you choose this option, unzip the file before proceeding to the next step.

### Compiling the C++ library
Here we show how to compile only the C++ part of Reaktoro. Its Python interface is an optional component of the project, and its compilation and installation is shown in the next section.

To build and install Reaktoro, ensure that [CMake](https://cmake.org/) is installed in your system. Reaktoro uses CMake for managing the whole build process, including the installation of third party libraries. 

Once CMake has been installed, go inside the directory of the downloaded Reaktoro source code. In the terminal, execute the following commands:
    
    mkdir build && cd build
    cmake ..
    make -j3

The first line above creates a directory called `build` and changes the current directory to it in Linux and MacOS systems. The command `cmake ..`  tells CMake to configure the build process based on the main `CMakeLists.txt` file in the root directory of Reaktoro's source code. Finally, `make -j3` compiles Reaktoro's source code using 3 parallel processes. To use all available processors in your machine, execute `make -j` instead. Be careful though, as the use of all available resources for compiling Reaktoro can freeze your machine!

To install the compiled libraries in your system, execute:
        
    make install

Note that this might require administrator rights, so that you would need to execute `sudo make install` instead. For a local installation, you can specify a directory path for the installed files using the CMake command:

    cmake .. -DCMAKE_INSTALL_PREFIX=/home/username/local/

### Compiling the Python interface
Most C++ classes and methods in Reaktoro are accessible from Python. To use its Python interface, Python wrappers to these C++ components must be compiled. These wrappers are generated using [Boost.Python](http://www.boost.org/doc/libs/1_60_0/libs/python/doc/html/index.html), so ensure your system has Boost installed, including `libboost_python`.

To build the Python wrappers, the CMake option `-DBUILD_PYTHON=ON` must be provided to the CMake command configuring the build process:

    cmake .. -DBUILD_PYTHON=ON

----------

*Compiling Reaktoro can take some time. This is because it heavily relies on template metaprogramming for efficient vector and matrix calculations, as well as for calculation of partial derivatives of most thermodynamic properties, such as activity coefficients, phase molar volumes, standard Gibbs energies, etc. In addition, it will also compile several third-party libraries, such as [CVODE](https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode) for efficient solution of ordinary differential equations (ODE), and the geochemical codes [PHREEQC](http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/) and [GEMS](http://gems.web.psi.ch/). Compilation of the Python wrappers can also take several minutes, as Boost.Python too relies on template metaprogramming.* 

----------

## Usage
We briefly show here some basic usage of Reaktoro for performing multiphase chemical equilibrium and kinetics calculations. More advanced use and customization is present in its user manual (work in progress!).

### Chemical equilibrium calculations

```c++
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database database("supcrt98.xml");

    ChemicalEditor editor(database);
    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- HCO3- CO2(aq) CO3--");
    editor.addGaseousPhase("H2O(g) CO2(g)");
    editor.addMineralPhase("Halite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(60, "celsius");
    problem.setPressure(300, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 100, "g");
    problem.add("NaCl", 1, "mol");

    ChemicalState state = equilibrate(problem);

    std::cout << state << std::endl;
}
```
### Chemical kinetics calculations

----------

## License

Reaktoro is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Reaktoro is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Reaktoro. If not, see <http://www.gnu.org/licenses/>.

----------

##Contact

For comments and requests, send an email to:

    allan.leal@erdw.ethz.ch
