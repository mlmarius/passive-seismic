```
$ python generate.py -h
usage: generate.py [-h] [--inv INVENTORY] [--cat CATALOG] [--outdata DATA]
                   [--outrf RFF]
                   asdffile

Calculation of receiver functions for all stations in a temporary survey..

positional arguments:
  asdffile         filepath of the asdf input file .

optional arguments:
  -h, --help       show this help message and exit
  --inv INVENTORY  filepath of the station xml file.
  --cat CATALOG    filepath of the event catalog xml file.
  --outdata DATA   filepath of the desired output waveform file in asdf
                   format.
  --outrf RFF      filepath of the desired output receiver function file in
                   asdf format.
```
