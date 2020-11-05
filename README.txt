To install all needed packages for _not_link, first install knuplot.
You find this in the Knuplot folder outside this folder.
Then:

pip install -r requirements_not_link.txt

(or sudo pip install -r requirements.txt if you were denied permission).

If you make symbolic links you will also need the following libraries, which
SINTEF has written, so you need them locally on your computer:

- Thermphys must be compiled, navigate to python-folder and install same way as Knuplot
- Thermopack must be compiled, navigate to pyThermopack (in addon folder) install with: python install/sudo python install
- libnum must me compiled, libsw depends on libnum
- libsw must be compiled, navigate to f2py_module (in addons folder), install
- by "python makescript.py optim" (I got an error messages here but it worked fine afterwards)
- trend must be compiled, go to pyTrend, build library: python makescript.py optim, install library: python install/sudo python install
