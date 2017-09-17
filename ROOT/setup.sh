#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/global/projecta/projectdirs/atlas/bnachman/code/pythia8226/
    export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}lib/:$LD_LIBRARY_PATH
}

setup_ROOT() {
    source /global/projecta/projectdirs/atlas/bnachman/code/root/bin/thisroot.sh
}

setup_fastjet() {
    export FASTJETLOCATION=/global/projecta/projectdirs/atlas/bnachman/code/fastjet-install/
    export LD_LIBRARY_PATH=/global/projecta/projectdirs/atlas/bnachman/code/fastjet-install/lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/global/projecta/projectdirs/atlas/bnachman/code/include/
    export BOOSTLIBLOCATION=/global/projecta/projectdirs/atlas/bnachman/code/lib/
    export LD_LIBRARY_PATH=${BOOSTLIBLOCATION}:$LD_LIBRARY_PATH
}

setup_ROOT
setup_PYTHIA
setup_fastjet
setup_boost

