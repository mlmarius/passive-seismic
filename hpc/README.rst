Creating a Python virtualenv on Raijin for passsive-seismic project
===================================================================

This is a quick guide to getting the passive-seismic software up and
running in both Portable Batch System (PBS) environment with MPI
support and virtual desktops (VDI). This setup is common in High Performance 
Compute (HPC) systems such as the National Computational Infrastructure's `Raijin
<http://nci.org.au/systems-services/national-facility/peak-system/raijin/>`__
system.

The instructions below are tailored to the NCI systems, including Raijin and 
virtual desktops (VDI). Some instructions may not be applicable depending on 
the setup of the HPC system you are using. They should apply to both single-node and
multi-node jobs; just set the number of cpus in the PBS directives in
the job submission script accordingly (e.g. ncpus=32 for 2 nodes).

These instructions assume you are using bash shell.

----------------
Pre-installation
----------------

Note: These instructions currently only work with gcc and not Intel compilers.

1. Clone the ``passive-seismic`` repository into your home directory, or
another directory
of your choice:

   .. code:: bash

       $ cd ~
       $ git clone https://github.com/GeoscienceAustralia/passive-seismic.git

2. Unload the icc compiler and default openmpi from the terminal:

   .. code:: bash

       $ module unload intel-cc
       $ module unload intel-fc
       $ module unload openmpi

3. Load the modules required for installation and running:

   **Python 3.5**
      .. code:: bash

       $ module load python3/3.4.3 python3/3.4.3-matplotlib
       $ module load hdf5/1.8.14p openmpi/1.8 mpi4py/2.0.0

   **Python 2.7**
      .. code:: bash

       $ module load python/2.7.3 python/2.7.3-matplotlib
       $ module load hdf5/1.8.14p openmpi/1.8.4 mpi4py/2.0.0-ompi-1.8.4-py2.7
       $ module load scipy/0.17.1-py2.7 cython/0.24-py2.7 matplotlib/1.5.1-py2.7 
       $ module load netcdf4-python/1.2.4-ncdf-4.3.3.1-py2.7 
       $ module load basemap/1.0.7-py2.7 pandas/0.18.1-py2.7

   (Alternatively, you may wish to add the above lines to your
   ``~/.profile`` file)

4. Now add the following lines to the end of your ``~/.profile`` file:

   **Python 3.5**
      .. code:: bash

       export PATH=$HOME/.local/bin:$PATH
       export PYTHONPATH=$HOME/.local/lib/python3.4/site-packages:$PYTHONPATH
       export VIRTUALENVWRAPPER_PYTHON=/apps/python3/3.4.3/bin/python3
       export LC_ALL=en_AU.UTF-8
       export LANG=en_AU.UTF-8
       source $HOME/.local/bin/virtualenvwrapper.sh

   **Python 2.7**
      .. code:: bash

       export PATH=$HOME/.local/bin:$PATH
       export PYTHONPATH=$HOME/.local/lib/python2.7.3/site-packages:$PYTHONPATH
       export VIRTUALENVWRAPPER_PYTHON=/apps/python/2.7.3/bin/python
       export LC_ALL=en_AU.UTF-8
       export LANG=en_AU.UTF-8
       source $HOME/.local/bin/virtualenvwrapper.sh

5. Install virtualenv and ``virtualenvwrapper`` in a terminal:

   **Python 3.5**
      .. code:: bash

       $ pip3 install  --user virtualenv virtualenvwrapper

   **Python 2.7**
      .. code:: bash

       $ pip install  --user virtualenv virtualenvwrapper

6. Refresh your environment by sourcing your ``~/.profile`` file:

   .. code:: bash

       $ source ~/.profile

------------
Installation
------------

1. Create a new virtualenv for ``passive-seismic``:

   .. code:: bash

       $ mkvirtualenv --system-site-packages seismic

2. Make sure the virtualenv is activated:

   .. code:: bash

       $ workon seismic

3. Clone ``h5py`` from ``https://github.com/basaks/h5py.git``:

   .. code:: bash

       $ cd ~
       $ git clone https://github.com/basaks/h5py.git
       $ cd ~/h5py
       $ export CC=mpicc
       $ python setup.py configure --mpi --hdf5=/apps/hdf5/1.8.14p
       $ python setup.py install


4. Install ``passive-seismic``:

   .. code:: bash

       $ cd ~/passive-seismic
       $ python setup.py install

5. Once installation has completed, you can run the tests to verify
   everything has gone correctly:

   .. code:: bash

       $ pip install pytest
       $ pytest tests/

