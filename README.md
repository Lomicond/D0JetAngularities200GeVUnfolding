# UnfoldingD0Jets

Analysis workflow for D0-meson-tagged jet unfolding.

## Environment

This analysis was tested in the STAR container environment with:

- ROOT 5.34.38
- gcc 4.8.5
- RooFit v3.60
- RooUnfold 2.0.1

Required RooUnfold version:

```bash
RooUnfold 2.0.1
commit d3526acf9540e1eac8e1e92cd26de21fa32c8ca6
repo https://github.com/roofit-dev/RooUnfold.git

```bash
setenv ROOUNFOLD /gpfs/mnt/gpfs01/star/pwg/lomicond/Ondrej/Jets/RooUnfold
setenv ROOT_INCLUDE_PATH ${ROOUNFOLD}/src:${ROOUNFOLD}:${ROOT_INCLUDE_PATH}
setenv LD_LIBRARY_PATH ${ROOUNFOLD}:${LD_LIBRARY_PATH}

```bash
cat > .rootlogon.C <<'EOF'
{
  gROOT->ProcessLine(".include /gpfs/mnt/gpfs01/star/pwg/lomicond/Ondrej/Jets/RooUnfold/src");
  gSystem->Load("/gpfs/mnt/gpfs01/star/pwg/lomicond/Ondrej/Jets/RooUnfold/libRooUnfold.so");
}
