if (${?FASTCALIBRATOR}) then
echo "already set"
else
setenv THISDIR `pwd`

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
# For CLHEP
setenv CLHEPSYS /afs/cern.ch/sw/lcg/external/clhep/2.0.4.2/slc4_amd64_gcc34/
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CLHEPSYS}/lib

if (${?DYLD_LIBRARY_PATH}) then
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/lib
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:${THISDIR}/../NtuplePackage/lib
endif



setenv PATH ${PATH}:${THISDIR}/bin
setenv PATH ${PATH}:${CLHEPSYS}/bin
setenv PATH ${PATH}:${THISDIR}/../NtuplePackage/bin

setenv NTUPLEPKGINCLUDE ${THISDIR}/../NtuplePackage/interface
setenv NTUPLEPKGLIB ${THISDIR}/../NtuplePackage/lib

setenv FASTCALIBRATOR ${THISDIR}/
endif
