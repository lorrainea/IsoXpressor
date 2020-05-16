#! /bin/sh

unzip isoSegmenter.zip

cd isoSegmenter
apt-get install libgd-dev libgif-dev
pip .
unzip REAL.zip
cd REAL
./configure
make
make install
mv ./src/real ../


