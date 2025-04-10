#!/bin/bash

# exit on any error
set -e

if [[ $# != 1 ]]; then
	echo "usage: $0 <new-repo-dir>"
	exit -1
fi

if [[ ! -f dimensions.yaml ]]; then
	echo "./dimensions.yaml file not found."
	echo "please run this script from the directory that contains it."
	exit -1
fi

REPO="$1"

if [[ -e "$REPO" ]]; then
	echo "$REPO already exists. cowardly refusing to proceed."
	exit -1
fi

##
## Now create the repo and the required dimensions
##

# Create repo
butler create "$REPO" --dimension-config=dimensions.yaml

# Patch butler.yaml
# This should be fixed in pipetasks (I think)
patch "$REPO/butler.yaml" butler.yaml.patch

# instruments we care about
butler register-instrument "$REPO" 'lsst.obs.lsst.LsstCam'
butler register-instrument "$REPO" 'lsst.obs.lsst.LsstComCam'
butler register-instrument "$REPO" 'lsst.obs.subaru.HyperSuprimeCam'

# regular datasets
butler register-dataset-type "$REPO" dia_source_visit    ArrowAstropy instrument visit

# SSP datasets
butler register-dataset-type "$REPO" sspDiaSourceInputs  ArrowAstropy instrument
butler register-dataset-type "$REPO" sspTrackletSources  ArrowAstropy instrument
butler register-dataset-type "$REPO" sspTrackletToSource ArrowAstropy instrument
butler register-dataset-type "$REPO" sspTracklets        ArrowAstropy instrument
butler register-dataset-type "$REPO" sspVisitInputs      ArrowAstropy instrument
butler register-dataset-type "$REPO" sspHypothesisTable  ArrowAstropy instrument
butler register-dataset-type "$REPO" sspEarthState       ArrowAstropy instrument
butler register-dataset-type "$REPO" sspLinkage          ArrowAstropy instrument sspHypothesisBundle
butler register-dataset-type "$REPO" sspLinkageSources   ArrowAstropy instrument sspHypothesisBundle

echo "Repository created in '$REPO'"
