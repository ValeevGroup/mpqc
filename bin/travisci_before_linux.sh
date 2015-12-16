#!/bin/sh

set -e

sudo add-apt-repository ppa:boost-latest/ppa -y
# Add PPA for a newer version GCC
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
# PPA for newer cmake (3.2.3)
sudo add-apt-repository ppa:george-edison55/precise-backports -y

sudo apt-get update -qq
sudo apt-get -y -qq --no-install-suggests --no-install-recommends --force-yes install cmake cmake-data
sudo apt-get install -y g++-$GCC_VERSION gfortran-$GCC_VERSION llvm-3.4 llvm-3.4-dev libgmp-dev libeigen3-dev libboost1.55-dev

cmake --version
gfortran-$GCC_VERSION --version
