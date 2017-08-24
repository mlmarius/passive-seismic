#!/usr/bin/env bash

if [ "$#"  -lt 1 ]; then
echo "Usage: $0 [mseed-volume]";
exit 0;
fi
set -x; seiscomp stop; seiscomp start spread scmaster
DBFLAG="mysql://sysop:sysop@localhost/seiscomp3"
VERBOSITY="-v"

scautopick --ep --playback -I file://$1 -d $DBFLAG > picks.xml
scautoloc --ep picks.xml -d $DBFLAG $VERBOSITY > origins.xml
scamp --ep origins.xml -I file://$1 -d $DBFLAG $VERBOSITY > amps.xml
scmag --ep amps.xml -d $DBFLAG $VERBOSITY > mags.xml
scevent --ep mags.xml -d $DBFLAG $VERBOSITY > events.xml

# scdb -i events.xml -d $DBFLAG $VERBOSITY
