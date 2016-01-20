# Reaktoro

Reaktoro is a unified framework for modeling chemically reactive systems. It provides methods for chemical equilibrium and kinetics calculations for multiphase systems. Reaktoro is mainly developed in C++ for performance reasons. A Python interface is available for a more convenient and simpler use of the scientific library. Currently Reaktoro can interface with two widely used geochemical software: [PHREEQC](http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/) and [GEMS](http://gems.web.psi.ch/). This document describes how to download and install Reaktoro, and demonstrate some basic usage.

## Installation

In the steps below we will show how one can download Reaktoro, build, and install it in Linux and MacOS systems. We plan to release binaries (i.e., the libraries already compiled) for Windows soon. Please get in touch so we can know how urgent these binaries are for you.

>**Note**: Compiling Reaktoro can take some time. This is because it heavily relies on template metaprogramming for efficient vector and matrix calculations, as well as for calculation of partial derivatives of most thermodynamic properties, such as activity coefficients, phase molar volumes, standard Gibbs energies, etc. In addition, it will also compile several third-party libraries, such as [CVODE](https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode) for efficient solution of ordinary differential equations (ODE), and the geochemical codes [PHREEQC](http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/) and [GEMS](http://gems.web.psi.ch/). Compilation of the Python wrappers can also take several minutes, as Boost.Python too relies on template metaprogramming.

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


## Usage
We briefly show here some basic usage of Reaktoro for performing multiphase chemical equilibrium and kinetics calculations. More advanced use and customization is present in its user manual (work in progress!).

### Chemical equilibrium calculations
Reaktoro contains methods for chemical equilibrium calculations based on either the *Gibbs energy minimization* (GEM) approach or on the *law of mass-action* (LMA) approach. In this section we describe step-by-step its use for performing chemical equilibrium calculations. 

#### A basic equilibrium calculation
In the code snippet below we show how the C++ interface of Reaktoro can be used to:

1. initialize a thermodynamic database;
2. specify the phases, and their species, composing the chemical system;
3. specify the equilibrium conditions for the calculation; and
4. perform the equilibrium calculation using a default method.

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
The first two lines:
```c++
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;
```
include the main Reaktoro header file: `Reaktoro.hpp`. By doing this, the application has access to all its classes and methods. The second line above is needed for convenience reasons: it eliminates the need to explicitly specify the namespace of Reaktoro components. Without it, we would need to write `Reaktoro::Database`, `Reaktoro::ChemicalSystem`,  `Reaktoro::EquilibriumProblem`, and so forth, which is a lot more verbose.

The equilibrium calculation uses the SUPCRT database together with the revised  *Helgeson-Kirkham-Flowers* (HKF) equations of state for the calculation of standard thermodynamic properties of aqueous, gaseous, and mineral species at temperatures 0 to 1000 °C and pressures 1 to 5000 bar. Thus, the following line is needed to initialize a `Database` instance using a database file `supcrt98.xml` that is found **at the same directory from where the application is executed!**  
```python
Database database("supcrt98.xml");
```

Once the `Database` instance has been initialized, one can use it to define the chemical system. For this, it is convenient to use the `ChemicalEditor` class, which currently permits the specification of aqueous, gaseous, and mineral phases, as well as specifying temperature and pressure interpolation points for the standard thermodynamic properties and choosing the equations of state (e.g., HKF, Pitzer, Peng-Robinson, and many others) for calculation of activity/fugacity coefficients of the species. In the lines below we use `ChemicalEditor`class to create an aqueous phase (with species separated by space!), a gaseous phase, and a single mineral phase. 
```python
ChemicalEditor editor(database);
editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ HCO3- CO2(aq) CO3--");
editor.addGaseousPhase("H2O(g) CO2(g)");
editor.addMineralPhase("Halite");
editor.addMineralPhase("Calcite");
```
The names of the species listed here must be found in the provided database file `supcrt98.xml`, otherwise an exception will be thrown. Note that there can be only one aqueous and one gaseous phase in the chemical system. Calling methods `addAqueousPhase` and `addGaseousPhase` more than once will simply replace the definition of those phases. However, method `addMineralPhase` can be called as many times as there are mineral phases in the system. If the mineral phase is a solid solution, then specify the mineral end-members like it was done in `addAqueousPhase` and `addGaseousPhase`. Note that a chemical system does not necessarily need to have an aqueous phase, or a gaseous phase, or mineral phases. Choose the combination of phases that describes your problem.

After the chemical system has been defined using class `ChemicalEditor`, it is now time to initialize an instance of `ChemicalSystem` class:

```python
ChemicalSystem system(editor);
``` 

The `ChemicalSystem` class is one of the most important classes in Reaktoro. It is the class used to computationally represent a chemical system, with all its phases, species, and elements. It is also the class used for computation of thermodynamic properties of phases and species, such as activities, chemical potentials, standard Gibbs energies, enthalpies, phase molar volumes, densities, and many others. Many classes in Reaktoro require an instance of `ChemicalSystem` for their initialization, since any chemical calculation needs to know the composition of the chemical system and the equations of state describing the non-ideal behavior of the phases.

Reaktoro provides the class `EquilibriumProblem` for convenient description of equilibrium conditions. Using this class allows one to set the temperature and pressure at equilibrium, and a recipe that describes a mixture of substances and their amounts, which can be seen as initial conditions for the equilibrium calculation:
```python
EquilibriumProblem problem(system);
problem.setTemperature(60, "celsius");
problem.setPressure(300, "bar");
problem.add("H2O", 1, "kg");
problem.add("CO2", 100, "g");
problem.add("NaCl", 1, "mol");
```
The units above can be changed, or even suppressed. If not provided, default units are used, such as K for temperatures, Pa for pressures, and mol for amounts. The `add` method in `EquilibriumProblem` supports both amount and mass units, such as `mmol`,  `umol`, `g`, `mg`, etc.

> **Note**: Make sure you use these units consistently. Using temperature units when setting pressure, for example, will throw an exception!

Once the equilibrium problem has been defined, it is now time to solve it. This can be done using the utility method `equilibrate`:
```c++
ChemicalState state = equilibrate(problem);
```
The line above uses the definition of the equilibrium problem stored in object `problem` to perform the equilibrium calculation. The result of the calculation is the object `state`, an instance of  `ChemicalState` class, which is used to store the chemical state (i.e., the temperature, pressure, and molar amounts of all species) of the system at prescribed equilibrium conditions. The `ChemicalState` class contains also methods for querying thermodynamic properties of the system.

Finally, we output the chemical state of the system to the standard output using:
```c++
std::cout << state << std::endl;
```
This will output tables describing the chemical state of the system. For example, the molar amounts, molar fractions, activities, activity coefficients, and chemical potentials of the species. The molar volumes of the phases, the amounts of each element in the phases, and also the total phase molar amounts.

### Chemical kinetics calculations

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


## Contact

For comments and requests, send an email to:

    allan.leal@erdw.ethz.ch