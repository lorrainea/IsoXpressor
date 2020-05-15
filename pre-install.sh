#! /bin/sh

unzip isoSegmenter.zip
unzip REAL.zip

cd REAL
./configure
make
make install
mv ./src/real ../


