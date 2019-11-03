Remote Debug of Reaktoro with PyDev
===================================

Make sure first that the PyDev plugin has been installed in Elipse.

The first step for remote debugging of a Reaktoro application is to inialize the `Debug Server` in PyDev. For this, go to `PyDev → Start Debug Server` on the top menu of Eclipse. You need to be in the *PyDev Perspective*, which can be open by following `Window → Open Perspective → Other` and choosing `PyDev`. If the option `Start Debug Server` can not be found, try to do this in the *Debug Perspective*.

Next, ensure the python module `pydevd.py` can be found by adding the directory where it exists to PYTHONPATH. For example, if eclipse is installed in `/opt`, then this path should be something like `/opt/eclipse/plugins/org.python.pydev_X/pysrc/`, where `X` contains the version of the PyDev plugin.

If executing from the terminal, do:

~~~
export PYTHONPATH = $PYTHONPATH:/opt/eclipse/plugins/org.python.pydev_X/pysrc/
~~~

If executing from eclipse, just go to `Window → Preferences → PyDev → Interpreters → Python Interpreters`, go to the tab `Libraries` and click on `New Folder`. Go to that directory above and click `OK`.
