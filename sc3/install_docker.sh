#!/usr/bin/env bash
DBFLAG="mysql://sysop:sysop@localhost/seiscomp3"
SEISCOMP_ROOT="/opt/seiscomp3"

sudo yum remove docker docker-common docker-selinux docker-engine
sudo yum install -y yum-utils
sudo yum-config-manager \
    --add-repo https://download.docker.com/linux/centos/docker-ce.repo
sudo yum makecache fast
sudo yum install -y docker-ce
sudo systemctl start docker

# clone and build docker image
git clone https://github.com/basaks/sc3-playback.git
cd ${HOME}/sc3-playback/tests/ && \
    sudo docker build -t sc3 .

# copy db from base machine
sqlite3 -batch -init $SEISCOMP_ROOT/share/db/sqlite3.sql \
    ${HOME}/sc3-playback/test.db .exit

# dump inventory and config xml
cd ${HOME}/sc3-playback && \
    scxmldump -fI -o inventory.xml -d ${DBFLAG} && \
    scxmldump -fC -o config.xml -d ${DBFLAG}

# Import the inventory
scdb --plugins dbsqlite3 \
    -d sqlite3://${HOME}/sc3-playback/test.db \
    -i ${HOME}/sc3-playback/inventory.xml

# Import the configuration
scdb --plugins dbsqlite3 \
    -d sqlite3://$HOME/sc3-playback/test.db \
    -i ${HOME}/sc3-playback/config.xml

# run the container as a daemon
sudo docker run -d --name sc3-ch -p 9999:22 \
    -v ${HOME}/sc3-playback/:/home/sysop/sc3-playback sc3

# Run using
# sudo docker exec -it --user sysop sc3-ch /bin/bash
