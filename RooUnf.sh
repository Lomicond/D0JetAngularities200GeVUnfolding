setup 64b && setup root 6.20.08
setenv ROOUNFOLD /gpfs/mnt/gpfs01/star/pwg/lomicond/Ondrej/Jets/RooUnfold
setenv ROOT_INCLUDE_PATH ${ROOUNFOLD}/src:${ROOUNFOLD}:${ROOT_INCLUDE_PATH}
setenv LD_LIBRARY_PATH ${ROOUNFOLD}:${ROOUNFOLD}/lib:${LD_LIBRARY_PATH}
