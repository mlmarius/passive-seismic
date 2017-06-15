Using the ANU Logfile Script
============================

The ``passive-seismic/convert_logs/decode_datfile.py`` can be used to convert
the binary `.dat` ANU log files into ``jsons``.

Can only use python2.7 for now.

------------------
Python Environment
------------------

Environment has these python packages:

   .. code:: bash

    click (6.7)
    numpy (1.13.0)
    pip (9.0.1)
    setuptools (36.0.1)
    wheel (0.29.0)

------------------
Running the script
------------------

There are two ways to run the ``anulog`` conversions:

#. directly run the python script
#. after ``passive-seismic`` is installed, you can use the command ``anulog``

Checkout the help string:

   .. code:: bash

    $ python decode_datfile.py --help

After ``passive-seismic`` is installed, you can simply use the ``anulog`` command

   .. code:: bash

    $ anulog --help

Both will produce the following help string:

::

  Usage: decode_datfile.py [OPTIONS] DATFILE
  Program to display contents of the logfile <datfile>.dat
  Options:
   -b, --bad_gps BOOLEAN      Print bad gps info
   -u, --gps_update BOOLEAN   Print the gps update info
   -i, --id_str BOOLEAN       Print bad id strings
   -t, --temperature BOOLEAN  Print bad temperature info
   -a, --all_print BOOLEAN    Print all
   -y, --year INTEGER RANGE   Gpsyear. max(Gpsyear - year) == 1
   -o, --output FILENAME      output json file name
   --help                     Show this message and exit.

A typical log file conversion command is just the following:

   .. code:: bash

    $ python decode_datfile.py logfile.dat -o output.json

This will output a ``output.json`` corresponding tot the ``logfile.dat``.

---------------------------
Using inside custom scripts
---------------------------

Once ``passive-seismic`` is installed, you have access to the decode
functionality ``decode_anulog`` to use in your script. Import it in your script
like the following:

   .. code:: bash

    $ In [1]: from convert_logs.decode_datfile import decode_anulog