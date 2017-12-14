#!/bin/sh

# delete all the annoying files,
#
find . -name '*.m~' -type f -delete
find . -name '*.sh~' -type f -delete
find . -name '*.txt~' -type f -delete
find . -name '*.tex~' -type f -delete
find . -name '.DS_Store' -type f -delete

cd beamforming
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..

cd china-passive
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..

cd dispersion-ftan
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..

cd dispersion-masw
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..

cd erb-passive
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..

cd gpr-benin
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..

cd interferometry
rm -rf data shared
ln -s ../data
ln -s ../shared
cd ..