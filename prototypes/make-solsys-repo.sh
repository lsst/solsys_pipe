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

if ! command -v butler >/dev/null 2>&1; then
	echo "butler is not on the path. Source and setup the LSST stack first."
	exit -1
fi

##
## Now create the repo and the required dimensions
##

# Create repo
butler create "$REPO" --dimension-config=dimensions.yaml

# Patch butler.yaml to add file access template for sspHypothesisBundle
# This should be changed upstream in:
#    https://github.com/lsst/daf_butler/blob/main/python/lsst/daf/butler/configs/datastores/fileDatastore.yaml#L8
# once we know what we want.
patch "$REPO/butler.yaml" butler.yaml.patch

# instruments we care about
butler register-instrument "$REPO" 'lsst.obs.lsst.LsstCam'
butler register-instrument "$REPO" 'lsst.obs.lsst.LsstComCam'
butler register-instrument "$REPO" 'lsst.obs.subaru.HyperSuprimeCam'

# regular datasets
butler register-dataset-type "$REPO" dia_source_detector ArrowAstropy instrument visit detector
butler register-dataset-type "$REPO" visit_summary ExposureCatalog band instrument day_obs physical_filter visit 

# SSP datasets
butler register-dataset-type "$REPO" sspEarthState       ArrowAstropy
butler register-dataset-type "$REPO" sspHypothesisTable  ArrowAstropy ssp_hypothesis_table
butler register-dataset-type "$REPO" sspDiaSourceInputs  ArrowAstropy instrument day_obs
butler register-dataset-type "$REPO" sspVisitInputs      ArrowAstropy instrument day_obs
butler register-dataset-type "$REPO" sspTrackletSources  ArrowAstropy instrument day_obs ssp_hypothesis_table
butler register-dataset-type "$REPO" sspTrackletToSource ArrowAstropy instrument day_obs ssp_hypothesis_table
butler register-dataset-type "$REPO" sspTracklets        ArrowAstropy instrument day_obs ssp_hypothesis_table
butler register-dataset-type "$REPO" sspLinkage          ArrowAstropy instrument day_obs ssp_hypothesis_table ssp_hypothesis_bundle
butler register-dataset-type "$REPO" sspLinkageSources   ArrowAstropy instrument day_obs ssp_hypothesis_table ssp_hypothesis_bundle
butler register-dataset-type "$REPO" sspBalancedLinkage  ArrowAstropy instrument day_obs ssp_hypothesis_table ssp_balanced_index
butler register-dataset-type "$REPO" sspBalancedLinkageSources ArrowAstropy instrument day_obs ssp_hypothesis_table ssp_balanced_index
echo "Repository created in '$REPO'"
