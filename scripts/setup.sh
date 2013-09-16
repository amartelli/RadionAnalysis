if [ -n "${FASTCALIBRATOR}" ] ; then
echo "already set"
else
export THISDIR=`pwd`

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THISDIR}/lib
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
# For CLHEP
export CLHEPSYS=/afs/cern.ch/sw/lcg/external/clhep/2.0.4.2/slc4_amd64_gcc34/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CLHEPSYS}/lib

if [ -n "${DYLD_LIBRARY_PATH}" ] ; then
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THISDIR}/lib
//export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
fi



export PATH=${PATH}:${THISDIR}/bin
export PATH=${PATH}:${CLHEPSYS}/bin
#export PATH=${PATH}:${THISDIR}/../NtuplePackage/bin

#export NTUPLEPKGINCLUDE=${THISDIR}/../NtuplePackage/interface
#export NTUPLEPKGLIB=${THISDIR}/../NtuplePackage/lib

export FASTCALIBRATOR=${THISDIR}/
fi
