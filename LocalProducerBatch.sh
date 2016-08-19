
# CERN batch script

cwd=$(pwd)

# ROOT 6
echo "Setting up environment"
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.18/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6/setup.sh

# Changing directory
cd /afs/cern.ch/user/m/mhaapale/work/public/LocalProducer

echo "Running ROOT script"
root mk_LocalProducer.C

echo "Success!"