# Using CMake to build a Reaktoro application in C++

This example demonstrates how your C++ application (executable) using Reaktoro
can be built with [CMake](https://cmake.org).

> Note: The description below has been tested using CMake v3.20.5. If you
> encounter errors, consider upgrading CMake or using alternative CMake
> configure/build commands.

In the root directory of this example, execute the following in the terminal:

~~~text
cmake -S . -B build
cmake --build build --config Release
~~~

If you built Reaktoro yourself, the commands above will fail because CMake
cannot find package Reaktoro. Assuming you built Reaktoro from the root
directory of the project with the commands:

~~~text
cd root/directory/of/reaktoro
cmake -P deps/install
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=installdir
cmake --build build --config Release --target install --parallel
~~~

then use the following commands instead:

~~~text
cmake -S . -B build -DCMAKE_PREFIX_PATH=root/directory/of/reaktoro/build/installdir
cmake --build build --config Release
~~~

> Note: Change `root/directory/of/reaktoro` to your own path (e.g.,
> `C:\Users\Username\codes\reaktoro` in Windows,
> `/home/username/codes/reaktoro` in Linux).

After building the application, you can run the application:

~~~text
./build/app
~~~

> Note: If you are using `conda` to manage the dependencies of Reaktoro when
> building it, then you should activate the Reaktoro conda environment
> **before** the CMake commands above:

~~~text
conda activate reaktoro
~~~

Read more on the installation instructions of Reaktoro.
