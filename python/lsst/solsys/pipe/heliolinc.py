import heliolinx.heliolinx as hl
import heliolinx.solarsyst_dyn_geo as sdg
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
from . import utils

#
# pipetask run -p ssp-heliolinc.yaml -b "$REPO" -i u/mjuric/test-small -o u/mjuric/test-small-output --register-dataset-types
# pipetask run -p ../pipelines/ssp-link.yaml -b "$REPO" -i u/mjuric/test-small -o u/mjuric/test-small-output -d "ssp_hypothesis_bundle.id in (1, 4)"
#
# This task searches tracklets for objects findable by a given bundle of hypothesis.
#
# The ssp_hypothesis_bundle ID is resolved to the actual hypotheses (the vector of numbers) by
# extracting the rows with the same ssp_hypothesis_bundle from a DataFrame stored as
# sspHypothesisDefinitions data type in the input collection.
#

__all__ = [
    "HeliolincTask",
    "HeliolincConfig",
]


class HeliolincConnections(lsst.pipe.base.PipelineTaskConnections,
                           dimensions=("instrument", "day_obs",
                                       "ssp_hypothesis_table", "ssp_hypothesis_bundle")):
    sspVisitInputs = connectionTypes.Input(
        doc="visit stats plus observer coordinates",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="visit_summary_dayobs_14"
    )
    sspTrackletSources = connectionTypes.Input(
        doc="sources that got included in tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_source_dayobs_14"
    )
    sspTracklets = connectionTypes.Input(
        doc="summary data for tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_dayobs_14"
    )
    sspTrackletToSource = connectionTypes.Input(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_to_source_dayobs_14"
    )
    sspHypothesisTable = connectionTypes.PrerequisiteInput(
        doc="hypotheses of asteroids' heliocentric radial motion",
        dimensions=["ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name = "sspHypothesisTable",
    )
    sspEarthState = connectionTypes.PrerequisiteInput(
        doc="Heliocentric Cartesian position and velocity for Earth",
        dimensions=[],
        storageClass="ArrowAstropy",
        name = "sspEarthState",
    )
    sspLinkage = connectionTypes.Output(
        doc="one line summary of each linkage",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_hypothesis_bundle", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_linkages",
    )
    sspLinkageSources = connectionTypes.Output(
        doc="indices connecting linkages (clusters) to trackletSources",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_hypothesis_bundle", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_linkage_sources",
    )

class HeliolincConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=HeliolincConnections):
    MJDref = lsst.pex.config.Field(
        dtype=float,
        default=0.0,
        doc="MJD or reference time. No sensible default is possible."
        )
    clustrad = lsst.pex.config.Field(
        dtype=float,
        default=1.0e5,
        doc="Clustering radius for the DBSCAN algorithm, in km."
        )
    clustchangerad = lsst.pex.config.Field(
        dtype=float,
        default=0.5,
        doc="Geocentric distance (AU), within which the clustering."
        )
    dbscan_npt = lsst.pex.config.Field(
        dtype=int,
        default=3,
        doc="Number of points npt for the DBSCAN algorithm"
        )
    minobsnights = lsst.pex.config.Field(
        dtype=int,
        default=3,
        doc="Minimum number of distinct observing nights for a valid linkage"
        )
    mintimespan = lsst.pex.config.Field(
        dtype=float,
        default=1.0,
        doc="Minimum timespan for a valid linkage, in days"
        )
    mingeodist = lsst.pex.config.Field(
        dtype=float,
        default=0.10,
        doc="Geocentric distance (AU) at the center of the innermost distance bin"
        )
    maxgeodist = lsst.pex.config.Field(
        dtype=float,
        default=100.0,
        doc="Minimum value in AU for the center of the outermost distance bin"
        )
    geologstep = lsst.pex.config.Field(
        dtype=float,
        default=1.5,
        doc="Factor by which distance increases from one bin to the next"
        )
    mingeoobs = lsst.pex.config.Field(
        dtype=float,
        default=0.0,
        doc="Minimum inferred geocentric distance for a valid tracklet"
        )
    minimpactpar = lsst.pex.config.Field(
        dtype=float,
        default=0.0,
        doc="Minimum inferred impact parameter (w.r.t Earth) for a valid tracklet"
        )
    use_univar = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Use the universal variable formulation of the Kepler equations, "
            "rather than the default fg function formulation"
        )
    max_v_inf = lsst.pex.config.Field(
        dtype=float,
        default=0.0,
        doc="Maximum value of v_infinity relative to the sun "
        "(must be greater than zero to probe interstellar orbits)"
        )
    autorun = lsst.pex.config.Field(
        dtype=int,
        default=1,
        doc="Automatically generate reference MJD instead of querying "
        "the user for it "
        )
    verbose = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Prints monitoring output."
        )

class HeliolincTask(lsst.pipe.base.PipelineTask):
    ConfigClass = HeliolincConfig
    _DefaultName = "link"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        bundle_id = butlerQC.quantum.dataId.ssp_hypothesis_bundle.id
        outputs = self.run(**inputs, ssp_hypothesis_bundle=bundle_id)
        butlerQC.put(outputs, outputRefs)

    def run(self, sspVisitInputs, sspTrackletSources, sspTracklets, sspTrackletToSource, sspHypothesisTable, sspEarthState, ssp_hypothesis_bundle):
        """doc string 
           here
        """

        # copy all config parameters from the Task's config object
        # to heliolinx's native config object.
        # Mask hypothesisTable to bundle_id == bundle_id

        config = hl.HeliolincConfig()
        allvars = [item for item in dir(hl.HeliolincConfig) if not item.startswith("_")]
        for var in allvars:
            setattr(config, var, getattr(self.config, var))

        sspHypothesisTable = sspHypothesisTable[sspHypothesisTable['bundle_id'] == ssp_hypothesis_bundle]
        sspHypothesisTable = sspHypothesisTable[['#r(AU)', 'rdot(AU/day)', 'mean_accel']]
        sspHypothesisTable.rename_columns(['#r(AU)', 'rdot(AU/day)', 'mean_accel'],
                                          ['HelioRad', 'R_dot', 'R_dubdot'])
        sspEarthState = sspEarthState[
            ['MJD', 'X', 'Y', 'Z', 'VX', 'VY', 'VZ']
        ]
        sspEarthState.rename_columns(['X', 'Y', 'Z', 'VX', 'VY', 'VZ'],
                                     ['x', 'y', 'z', 'vx', 'vy', 'vz'])

        (sspLinkage, sspLinkageSources) = hl.heliolinc(config,
                                                       utils.df2numpy(sspVisitInputs,      "hlimage"),
                                                       utils.df2numpy(sspTrackletSources,  "hldet"),
                                                       utils.df2numpy(sspTracklets,        "tracklet"),
                                                       utils.df2numpy(sspTrackletToSource, "longpair"),
                                                       utils.df2numpy(sspHypothesisTable,  "hlradhyp"),
                                                       utils.df2numpy(sspEarthState,       "EarthState")
                                                      )

        return lsst.pipe.base.Struct(sspLinkage=sspLinkage,
                                     sspLinkageSources=sspLinkageSources
                                     )
