#! /bin/sh

unzip isoSegmenter.zip

cd isoSegmenter
pip .
unzip REAL.zip
cd REAL
./configure
make
make install
mv ./src/real ../


