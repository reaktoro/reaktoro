
Follow the instructions [here](https://msys2.github.io/) to install MSYS2 and update its system packages.

1. Install 32-bit Python 2.7.x into C:\Python27-w32 without updating the PATH variable
2. Install 64-bit Python 2.7.x into C:\Python27-w64 without updating the PATH variable
3. Execute `C:\Python27-w32\Scripts\easy_install.exe -U pyinstaller pypiwin32`
4. Execute `C:\Python27-w64\Scripts\easy_install.exe -U pyinstaller pypiwin32`

C:\Python27-w64\Scripts\pip.exe install https://github.com/pyinstaller/pyinstaller/archive/develop.zip
C:\Python27-w32\Scripts\pip.exe install https://github.com/pyinstaller/pyinstaller/archive/develop.zip

Looks like Nuitka only is installed correctly with pip, not easy_install - the later results in some ImportError of pkg_resources.


The above will be used to create an executable for Reaktoro application that does not require installation of Python.

~~~
pacman -S cmake patch

pacman -S mingw-w64-{i686,x86_64}-toolchain

pacman -S mingw64/mingw-w64-x86_64-boost
pacman -S mingw64/mingw-w64-x86_64-python2
pacman -S mingw64/mingw-w64-x86_64-python3
pacman -S mingw64/mingw-w64-x86_64-python2-numpy
pacman -S mingw64/mingw-w64-x86_64-python2-pip
pacman -S mingw64/mingw-w64-x86_64-python2-pywin32

pacman -S mingw32/mingw-w64-i686-boost
pacman -S mingw32/mingw-w64-i686-python2
pacman -S mingw32/mingw-w64-i686-python3
pacman -S mingw32/mingw-w64-i686-python2-numpy
pacman -S mingw32/mingw-w64-i686-python2-pip
pacman -S mingw32/mingw-w64-i686-python2-pywin32

pacman -S rsync

cd C:\msys64\home\leal_a\Downloads\MINGW-packages\mingw-w64-python-pywin32

~~~

Install GnuWin32
================

Download [GnuWin32](http://sourceforge.net/projects/getgnuwin32/files/latest/download?source=files), execute it and extract it to a directory of your choice. Open the terminal, navigate to the directory where GnuWin32 was extracted and execute `download.bat`. This will download all GnuWin32 packages, which can take a long time. Once this is finished, run `install C:\gnuwin32`, which will install all downloaded packages to `C:\gnuwin32`. Add the directory `C:\gnuwin32\bin` to the system path by editing the environmental variable `Path` so that the GNU applications can be easily accessed from the terminal.

To avoid any potential issue, restart Windows.

Building 64-bit Reaktoro
========================

Open the Windows terminal and change directory to where Reaktoro was downloaded.

set PATHBKP=%PATH%

set PATH=C:\msys64\bin;C:\msys64\mingw64\bin;%PATH%

cmake -DCMAKE_PREFIX_PATH=C:\msys64\mingw64 -DPYTHON_LIBRARY=C:\msys64\mingw64\lib\libpython2.7.dll.a -DPYTHON_INCLUDE_DIR=C:\msys64\mingw64\include\python2.7 -DBUILD_ALL=ON -G"MSYS Makefiles" ../..





Make sure
Start -> MSYS2 64bit -> MinGW-w64 Win32 Shell

Make sure the mingw-w64 Python executable is accessible, not the one downloaded from the officeial Python web-page.

#-------------------------------------------------------------------------------
Execute C:\msys64\mingw32_shell.bat


cmake -DCMAKE_PREFIX_PATH=/mingw32/ -DPYTHON_LIBRARY=/mingw32/lib/libpython2.7.dll.a -DPYTHON_INCLUDE_DIR=/mingw32/include/python2.7 -DBUILD_ALL=ON -G"MSYS Makefiles" ../..

#-------------------------------------------------------------------------------
Execute C:\msys64\mingw64_shell.bat


cmake -DCMAKE_PREFIX_PATH=/mingw64/ -DPYTHON_LIBRARY=/mingw64/lib/libpython2.7.dll.a -DPYTHON_INCLUDE_DIR=/mingw64/include/python2.7 -DBUILD_ALL=ON -G"MSYS Makefiles" ../..



C:\Python27\;C:\Python27\Scripts;%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem;%SYSTEMROOT%\System32\WindowsPowerShell\v1.0\;C:\Program Files\PSI Tools\;%systemroot%\System32\WindowsPowerShell\v1.0\;C:\Program Files\gnuplot\bin;C:\Program Files (x86)\Windows Kits\8.1\Windows Performance Toolkit\;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files (x86)\CMake\bin;C:\Program Files (x86)\Git\cmd;C:\Program Files (x86)\Skype\Phone\


C:\Program Files\PSI Tools\;%systemroot%\System32\WindowsPowerShell\v1.0\;C:\Program Files\gnuplot\bin;C:\Program Files (x86)\Windows Kits\8.1\Windows Performance Toolkit\;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files (x86)\Skype\Phone\



C:\msys64\mingw64\bin;%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem;%SYSTEMROOT%\System32\WindowsPowerShell\v1.0\;C:\Program Files\PSI Tools\;%systemroot%\System32\WindowsPowerShell\v1.0\;C:\Program Files\gnuplot\bin;C:\Program Files (x86)\Windows Kits\8.1\Windows Performance Toolkit\;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files (x86)\CMake\bin;C:\Program Files (x86)\Git\cmd;C:\Program Files (x86)\Skype\Phone\




cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/reaktoro/release/ -DCMAKE_LIBRARY_PATH=/mingw64/lib -DBUILD_GEMS=ON -DBUILD_PHREEQC=ON -DBUILD_PYTHON_WRAPPERS=ON -G"MSYS Makefiles" ../..



pacman -Ss | grep '^mingw64/' -A 1 | sed -r -e 's/\[installed.*\]//' -e 's#^mingw64/([^ ]+)#<br/>mingw/<b>\1</b>#' -e 's/-x86_64//'

pacman -Ss | grep '^mingw64/' -A 1 | sed -r -e 's/\[installed.*\]//' -e 's#^mingw64/([^ ]+)#mingw/\1#' -e 's/-xi686//'
