Trouble shooting
==================

``OpenFabrics device`` error
------------------------------

When using the mpirun, you may encounter with::

   --------------------------------------------------------------------------
   By default, for Open MPI 4.0 and later, infiniband ports on a device
   are not used by default.  The intent is to use UCX for these devices.
   You can override this policy by setting the btl_openib_allow_ib MCA parameter
   to true.

   Local host:              lxbk1039
   Local adapter:           mlx5_0
   Local port:              1

   --------------------------------------------------------------------------
   --------------------------------------------------------------------------
   WARNING: There was an error initializing an OpenFabrics device.

   Local host:   lxbk1039
   Local device: mlx5_0
   --------------------------------------------------------------------------
   [lxbk1039:106838] 11 more processes have sent help message help-mpi-btl-openib.txt / ib port not selected
   [lxbk1039:106838] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
   [lxbk1039:106838] 11 more processes have sent help message help-mpi-btl-openib.txt / error in device init

You can add ``export OMPI_MCA_btl=^openib`` to your slurm file or system environment.
Please refer to `There was an error initializing an OpenFabrics device. #10693 <https://github.com/open-mpi/ompi/issues/10693>`_ for more information.

``ibv_exp_query_device`` error
-------------------------------

When using the mpirun, you may encounter with::

    ibv_exp_query_device: invalid comp_mask !!! (comp_mask = 0x7f81fdd30040 valid_mask = 0x1)
    [r4-node9][[44333,1],23][btl_openib_component.c:1670:init_one_device] 
    error obtaining device attributes for mlx5_0 errno says Invalid argument
    --------------------------------------------------------------------------
    WARNING: There was an error initializing an OpenFabrics device.
    Local host:   r4-node9
    Local device: mlx5_0
    --------------------------------------------------------------------------

You can run with ``mpirun --mca pml ucx ...`` to force the UCX PML to be used in MPI.
Please refer to this `ibv_exp_query_device error when runing mpi program using roce v2 #5807 <https://github.com/open-mpi/ompi/issues/5807>`_ for more information.


Notes on gcc Installation
-------------------------

Installation of gcc can be easily through HomeBrew in MacOS::

    brew install gcc
    
or through apt in Linux (Fedora/RHEL-based distros)::
    
    apt install build-essential

or through pacman in Windows::

    pacman -Syu
    pacman -S --needed base-devel mingw-w64-x86_64-toolchain

But sometimes you need to install gcc manually (in a HPC), the following steps may help you:

Before install gcc, you need to install the reliable lib: `gmp` , `mpfr`, `mpc`, and make sure that your exist environment unset, if you used set a environment in ``~/.bashrc``, ``~/.zshrc``,  or ``~/.bash_profile`` like::

    export CPLUS_INCLUDE_PATH={path to cpp include}:$CPLUS_INCLUDE_PATH

Then, you need to run ``unset CPLUS_INCLUDE_PATH`` in your terminal.

Now let's create a direction to store the local library firstly::

    mkdir $PATH_TO_LIBRARY

Replace ``$PATH_TO_LIBRARY`` with the path to where you want to install the libaray locally.

Installing gmp
>>>>>>>>>>>>>>

First, download the GMP library source code. You can get it from the official website or by using the `wget` command::

    wget https://gmplib.org/download/gmp/gmp-6.1.0.tar.bz2 --no-check-certificate


Extract the downloaded file::

    tar -xjf gmp-6.1.0.tar.bz2

Change to the extracted directory::

    cd gmp-6.1.0


Now, you need to decide on the installation directory. Let's say you want to install GMP in ``$PATH_TO_LIBRARY``. You can set the ``--prefix`` option to specify the installation directory. Run the following commands::

    ./configure --prefix=$PATH_TO_LIBRARY
    make

After the compilation process is complete, run the following command to check if the library has been compiled without any errors::

    make check

If everything went well, you can proceed with the installation::

    sudo make install

Installing mpfr
>>>>>>>>>>>>>>>

Similarly, download the MPFR library source code::

   wget https://www.mpfr.org/mpfr-4.1.0/mpfr-4.1.0.tar.gz --no-check-certificate
   
Extract the downloaded file::

   tar -xzf mpfr-4.1.0.tar.gz

Change to the extracted directory::

   cd mpfr-4.1.0

Now, you need to decide on the installation directory. Let's say you want to install GMP in ``$PATH_TO_LIBRARY``. You can set the ``--prefix`` option to specify the installation directory. Run the following commands::

   
   ./configure --prefix=$PATH_TO_LIBRARY --with-gmp=$PATH_TO_LIBRARY
   make
   
After the compilation process is complete, run the following command to check if the library has been compiled without any errors::
   
   make check
   
If everything went well, you can proceed with the installation::

   sudo make install

Installing mpc
>>>>>>>>>>>>>>

Similarly, download the MPC library source code::

   wget https://ftp.gnu.org/gnu/mpc/mpc-1.2.1.tar.gz --no-check-certificate

Extract the downloaded file::

   tar -xzf mpc-1.2.1.tar.gz

Change to the extracted directory::

   cd mpc-1.2.1

Now, you need to decide on the installation directory. Let's say you want to install GMP in ``$PATH_TO_LIBRARY``. You can set the ``--prefix`` option to specify the installation directory. Run the following commands::

   ./configure --prefix=$PATH_TO_LIBRARY --with-mpfr=$PATH_TO_LIBRARY --with-gmp=$PATH_TO_LIBRARY
   make

After the compilation process is complete, run the following command to check if the library has been compiled without any errors::

   make check
   
If everything went well, you can proceed with the installation::

   sudo make install
   
Installing gcc
>>>>>>>>>>>>>>>>>>>>>>

Similarly, download the gcc library source code::

   wget https://ftp.gnu.org/gnu/gcc/gcc-11.2.0/gcc-11.2.0.tar.xz --no-check-certificate
   
Making direction to build gcc library::

   mkdir gcc-11.2.0-build
   cd gcc-11.2.0-build/
   
Now, let's configure the gcc to the ``$PATH_TO_LIBRARY``::

   ../gcc-11.2.0/configure --prefix=$PATH_TO_LIBRARY --disable-multilib --enable-languages=c,c++,fortran --with-mpc=$PATH_TO_LIBRARY --with-mpfr=$PATH_TO_LIBRARY --with-gmp=$PATH_TO_LIBRARY
   make -j 8
   
Install it::

   make install
   
To use the installed gcc library, you might need to add the installation path to your ``$LD_LIBRARY_PATH`` and ``$PATH`` environment variables::

   export PATH=/data01/home/zhangbo/library/bin:$PATH
   export LD_LIBRARY_PATH=/data01/home/zhangbo/library/lib64:$LD_LIBRARY_PATH

To make these changes permanent, add the above lines to your shell's configuration file (e.g., ``~/.bashrc``, ``~/.zshrc``, or ``~/.bash_profile``).
