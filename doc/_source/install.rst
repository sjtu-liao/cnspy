Getting Started
===============

Installation
------------

From Github Repository::

    cd $INSTALL_LOCATION
    git clone https://github.com/Zilpher/cnspy.git

Replace ``$INSTALL_LOCATION`` with the path to where you want to install the code.

Using `conda` to create a virture environment and install the necessery package::

    conda create -n <envname> numpy sympy scipy jupyter matplotlib scipy

Then, activate this environment::

    conda activate <envname>

To install the code, run::

    python setup.py bdist_wheel
    pip install dist/cnspy-0.1-py3-none-any.bdist_wheel    

To uninstall, run::

    pip uninstall cnspy

To reinstall, run::

    pip install dist/cnspy-0.1-py3-none-any.whl --force-reinstall

Prerequisites\:

* Python 3.9 
* Python packages

    * NumPy
    * SciPy
    * Sympy
    * Matplotlib

* gcc 4.8.5

    * gmp 6.2.0
    * mpfr 4.2.0

* openmpi 2.1.0
