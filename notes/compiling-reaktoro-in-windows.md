
Follow the instructions [here](https://msys2.github.io/) to install MSYS2 and update its system packages. 

~~~
pacman -S cmake patch

pacman -S mingw-w64-{i686,x86_64}-toolchain
pacman -S mingw64/mingw-w64-{i686,x86_64}-boost
pacman -S mingw64/mingw-w64-{i686,x86_64}-python2
pacman -S mingw64/mingw-w64-{i686,x86_64}-python3
pacman -S mingw64/mingw-w64-{i686,x86_64}-python2-numpy
~~~

Make sure 
Start -> MSYS2 64bit -> MinGW-w64 Win32 Shell

Make sure the mingw-w64 Python executable is accessible, not the one downloaded from the officeial Python web-page.

#-------------------------------------------------------------------------------
Execute C:\msys64\mingw32_shell.bat

cmake -DCMAKE_INSTALL_PREFIX=/mingw32/ -DCMAKE_PREFIX_PATH=/mingw32/ -DPYTHON_LIBRARY=/mingw32/lib/libpython2.7.dll.a -DPYTHON_INCLUDE_DIR=/mingw32/include/python2.7 -DBUILD_GEMS=ON -DBUILD_PHREEQC=ON -DBUILD_PYTHON_WRAPPERS=ON -G"MSYS Makefiles" ../..

-DPYTHON_INCLUDE_DIR=C:/msys64/mingw32/include/python2.7
-DPYTHON_LIBRARY=C:/msys64/mingw32/lib/python2.7/config
#-------------------------------------------------------------------------------
Execute C:\msys64\mingw64_shell.bat

cmake -DCMAKE_INSTALL_PREFIX=/mingw64/ -DCMAKE_PREFIX_PATH=/mingw64/ -DBUILD_GEMS=ON -DBUILD_PHREEQC=ON -DBUILD_PYTHON_WRAPPERS=ON -G"MSYS Makefiles" ../..



C:\msys64\mingw64\bin;C:\Python27\;C:\Python27\Scripts;%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem;%SYSTEMROOT%\System32\WindowsPowerShell\v1.0\;C:\Program Files\PSI Tools\;%systemroot%\System32\WindowsPowerShell\v1.0\;C:\Program Files\gnuplot\bin;C:\Program Files (x86)\Windows Kits\8.1\Windows Performance Toolkit\;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files (x86)\CMake\bin;C:\Program Files (x86)\Git\cmd;C:\Python27;C:\Python27\Scripts;C:\Program Files\mingw-w64\x86_64-5.1.0-posix-seh-rt_v4-rev0\mingw64\bin;C:\Program Files (x86)\Skype\Phone\


C:\Program Files\PSI Tools\;%systemroot%\System32\WindowsPowerShell\v1.0\;C:\Program Files\gnuplot\bin;C:\Program Files (x86)\Windows Kits\8.1\Windows Performance Toolkit\;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files (x86)\Skype\Phone\



C:\msys64\mingw64\bin;%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem;%SYSTEMROOT%\System32\WindowsPowerShell\v1.0\;C:\Program Files\PSI Tools\;%systemroot%\System32\WindowsPowerShell\v1.0\;C:\Program Files\gnuplot\bin;C:\Program Files (x86)\Windows Kits\8.1\Windows Performance Toolkit\;C:\Program Files\Microsoft SQL Server\110\Tools\Binn\;C:\Program Files (x86)\Microsoft SDKs\TypeScript\1.0\;C:\Program Files (x86)\CMake\bin;C:\Program Files (x86)\Git\cmd;C:\Program Files (x86)\Skype\Phone\




cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/reaktoro/release/ -DCMAKE_LIBRARY_PATH=/mingw64/lib -DBUILD_GEMS=ON -DBUILD_PHREEQC=ON -DBUILD_PYTHON_WRAPPERS=ON -G"MSYS Makefiles" ../..