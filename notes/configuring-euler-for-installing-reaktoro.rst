Configuring Euler for installation of Reaktoro
==============================================

Installing Python
-----------------

~~~
cd ~/downloads
wget https://www.python.org/ftp/python/2.7.10/Python-2.7.10.tar.xz
tar -xvf Python-2.7.10.tar.xz
cd Python-2.7.10
./configure --prefix=$HOME/local --enable-shared
make -j 32
make install
~~~

Installing Python Packages
--------------------------

First, `pip` needs to be installed.

~~~
cd ~/downloads
wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py
~~~

Now that `pip` has been installed, any Python package can be installed via `pip install package-name`. Reaktoro needs the following packages:

~~~
pip install numpy pyinstaller tabulate
~~~

Installing CMake
----------------

The available CMake in Euler cluster does not support https protocol, which is needed when dowloading Reaktoro's third-party libraries using `ExternalProject_Add` command. To solve this, a local CMake needs to be built with support to https protocol. This can be done as:

~~~
cd ~/downloads
wget http://www.cmake.org/files/v3.3/cmake-3.3.1.tar.gz
tar -xvf cmake-3.3.1.tar.gz
cd cmake-3.3.1
./bootstrap
cmake . -DCMAKE_USE_OPENSSL=ON
make -j 32
make install
~~~

Installing Boost
----------------

~~~
cd ~/downloads
wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz
tar -xvf boost_1_59_0.tar.gz
cd boost_1_59_0
./bootstrap.sh --with-libraries=python --prefix=$HOME/local
./b2 -j 32 install
~~~

cmake ../.. -DCMAKE_PREFIX_PATH=~/local -DPYTHON_LIBRARY=/cluster/apps/python/2.7.6/x86_64/lib64/libpython2.7.so


module load gcc

~~~
module load python/2.7
pip install --user pyaml tabulate
~~~

The `--user` parameter will make `pip` to install the Python packages locally at `~/.local`.

Setting the .bashrc file
------------------------

~~~
# Set the paths to the directories where local binaries have been installer
export PATH=$HOME/local/bin/:$PATH

# Set the paths to the directories where local Python packages have been installed
export PATH=$HOME/.local/bin/:$PATH

# Load the Euler modules needed to build Reaktoro
module load gcc/4.9.2
module load python/2.7
module load boost/1.57.0
~~~
