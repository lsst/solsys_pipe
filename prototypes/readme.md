# Notes & cheat-sheets while developing

Assuming the stack will be installed in `~/projects/lsst`

## Developing pipelines at the USDF:
See https://rsp.lsst.io/v/usdfprod/guides/notebooks/science-pipelines/science-pipelines-development.html

## Installing the stack:

Follow the install procedure from:

	https://pipelines.lsst.io/install/lsstinstall.html

Note: use `./lsstinstall -T w_latest`, not the v_... tag

Then do:

```
eups distrib install -t w_latest lsst_distrib
```

setup lsst_distrib
conda activate lsst-scipipe-8.0.0
mamba install -c conda-forge ipython
mamba install -c conda-forge ipykernel

## Activating the environment

```
ROOT=~/projects/lsst
cd "$ROOT"
source loadLSST.bash
setup lsst_distrib
```

## Install solsys_pipe
```
PIPESRC=~/projects/github.com/lsst/solsys_pipe
# rm -rf solsys_pipe
git clone https://github.com/lsst/solsys_pipe
pushd solsys_pipe
git checkout u/mjuric/pipeline
setup -r .
scons
```

## Creating a new Butler repo

> From https://doc.lsst.eu/tutorial/butler.html
> To create new dimensions, modify this:
>    https://github.com/lsst/daf_butler/blob/main/python/lsst/daf/butler/configs/dimensions.yaml
> This is called "the universe" of dimensions
> More demos: https://github.com/LSSTScienceCollaborations/StackClub/blob/master/Basics/Gen3ButlerTutorial.ipynb
>

```
REPO=$ROOT/repo
pushd $PIPESRC/prototypes
butler create "$REPO" --dimension-config=dimensions.yaml
popd

butler register-instrument "$REPO" 'lsst.obs.lsst.LsstCam'
butler register-instrument "$REPO" 'lsst.obs.subaru.HyperSuprimeCam'

butler register-dataset-type "$REPO" sspDiaSourceInputs  DataFrame instrument
butler register-dataset-type "$REPO" sspTrackletSources  DataFrame instrument
butler register-dataset-type "$REPO" sspTrackletToSource DataFrame instrument
butler register-dataset-type "$REPO" sspTracklets        DataFrame instrument
butler register-dataset-type "$REPO" sspVisitInputs      DataFrame instrument
butler register-dataset-type "$REPO" sspHypothesisTable  DataFrame instrument
butler register-dataset-type "$REPO" sspEarthState       DataFrame instrument
butler register-dataset-type "$REPO" sspLinkage          DataFrame instrument sspHypothesisBundle
butler register-dataset-type "$REPO" sspLinkageSources   DataFrame instrument sspHypothesisBundle
butler query-dataset-types "$REPO" "ssp*" -v

pushd $PIPESRC/prototypes/data
butler ingest-files "$REPO" sspDiaSourceInputs u/mjuric/test-small sspDiaSourceInputs.csv
butler ingest-files "$REPO" sspVisitInputs     u/mjuric/test-small sspVisitInputs.csv
butler ingest-files "$REPO" sspHypothesisTable u/mjuric/test-small sspHypothesisTable.csv
butler ingest-files "$REPO" sspEarthState      u/mjuric/test-small sspEarthState.csv
popd

butler query-datasets "$REPO" "ssp*"

# Run the code from task-demo to inject some dimensions entries
# TODO: turn this into an executable!

butler query-dimension-records "$REPO" sspHypothesisBundle

### Build heliolinc
cd ~/projects/github.com/lsst-dm/heliolinc2
git checkout u/mjuric/pipeline
setup -r .
cd python
c++ -O3 -Wall -shared -fopenmp -std=c++14 -fPIC $(python3 -m pybind11 --includes) heliohypy.cpp ../src/solarsyst_dyn_geo01.cpp -o heliohypy$(python3-config --extension-suffix)
cd ../notebooks
pipetask run -p ssp-heliolinc.yaml -b "$REPO" -i u/mjuric/test-small -o u/mjuric/test-small-output --register-dataset-types

### With new dimensions
REPO=$HOME/projects/lsst_ssp/repo2
cd ~/projects/github.com/lsst-dm/heliolinc2/notebooks
butler create "$REPO" --dimension-config=dimensions.yaml

### solsys_pipe (kernel setup)
python -m ipykernel install --user --name mjuric_lsst --display-name "LSST (solsys_pipe)"
edit:
/astro/users/mjuric/.local/share/jupyter/kernels/mjuric_lsst/kernel.sh
/astro/users/mjuric/.local/share/jupyter/kernels/mjuric_lsst/lsst_kernel.sh

### Setup
cd ~/projects/lsst_ssp
conda activate lsst-scipipe-8.0.0
setup lsst_distrib
REPO=$HOME/projects/lsst_ssp/repo

## HL Ari chat

cols = (0,1,2)
accelmat = np.loadtxt('heliohyp_ombtest.txt', usecols=cols, skiprows=1, dtype=solardg.hlradhyp)
