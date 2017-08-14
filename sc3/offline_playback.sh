#!/usr/bin/env bash
if [ "$#"  -lt 2 ]; then
    echo "Usage: $0 [mseed-volume] [output-xml]";
    exit 0;
fi

set -x; seiscomp stop; seiscomp start spread scamaster

DBFLAG="mysql://sysop:sysop@localhost/seiscomp3"
CONFIGFLAGS="--verbosity=4"

FLAGS="$CONFIGFLAGS $DBFLAG"

echo "Cleaning Database"; seiscomp exec scdbstrip $FLAGS --days 0
echo "Starting autoloc..."; seiscomp exec scautoloc $FLAGS --playback --start-stop-msg=1 --auto-shutdown=1 --shutdown-master-module=scautopick &
echo "Starting amplitude and magtool..."

seiscomp exec scamp $FLAGS --start-stop-msg=1 --auto-shutdown=1 --shutdown-master-module=scautoloc &
seiscomp exec scmag $FLAGS --start-stop-msg=1 --auto-shutdown=1 --shutdown-master-module=scamp &
echo "Starting eventtool..."; seiscomp exec scevent --db-disable $FLAGS --start-stop-msg=1 --auto-shutdown=1 --shutdown-master-module=scmag &

echo "Starting sceplog..."; seiscomp exec sceplog $CONFIGFLAGS --auto-shutdown=1 --shutdown-master-module=scevent > $2 &
pid=$!
echo "Starting autopick..."; seiscomp exec scautopick --playback -I $1 $FLAGS --start-stop-msg=1

echo "Finished waveform processing - wait for finishing event processing"; wait $pid
echo "Finished event processing"
