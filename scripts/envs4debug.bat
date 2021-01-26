@REM This file will be copied by cmake to the root of the build directory.
@REM Execute `envs4debug` from that directory and the environment variables
@REM below will be update so that the python package and C++ libraries can
@REM be found if necessary (e.g., to execute pytest).

set PYTHONPATH=%CD%\python\package\installed\Debug\Lib\site-packages;%CD%\..\deps\build\install\public\lib\site-packages;%PYTHONPATH%
set PATH=%CD%\bin\Debug;%CD%\..\deps\build\install\public\lib;%PATH%
