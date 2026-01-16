import heliolinx.heliolinx as hl
import heliolinx.solarsyst_dyn_geo as sdg
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes, NoWorkFound
from . import utils

# This task purifies and de-duplicates candidate asteroid linkages output
# by heliolinc

class LinkPurifyConnections(lsst.pipe.base.PipelineTaskConnections,
                            dimensions=("instrument", "day_obs",
                                        "ssp_hypothesis_table", "ssp_balanced_index")):
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
        name="ssp_tracklet_sources"
    )
    sspBalancedLinkages = connectionTypes.Input(
        doc="one line summary of each linkage",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_balanced_index", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_balanced_linkages",
    )
    sspBalancedLinkageSources = connectionTypes.Input(
        doc="indices connecting linkages (clusters) to trackletSources",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_balanced_index", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_balanced_linkage_sources",
    )
    sspPurifiedLinkages = connectionTypes.Output(
        doc="one line summary of each purified linkage",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_balanced_index", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_purified_linkages",
    )
    sspPurifiedLinkageSources = connectionTypes.Output(
        doc="indices connecting linkages (clusters) to trackletSources",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_balanced_index", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_purified_linkage_sources",
    )

class LinkPurifyConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=LinkPurifyConnections):
    useorbMJD = lsst.pex.config.Field(
        dtype=int,
        default=1,
        doc="Use the MJD of the orbit rather than the reference MJD, if available"
    )
    simptype = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Type of simplex used to initialize orbit fitting"
    )
    ptpow = lsst.pex.config.Field(
        dtype=int,
        default=-1,
        doc="Exponent used for number of points in linkage quality metric."
            + "Negative activates product of per-night point counts."
    )
    nightpow = lsst.pex.config.Field(
        dtype=int,
        default=-1,
        doc="Exponent used for number of nights in linkage quality metric."
            + "Negative activates product of per-night point counts."
    )
    timepow = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Exponent used for linkage time span in linkage quality metric"
    )
    rmspow = lsst.pex.config.Field(
        dtype=int,
        default=1,
        doc="Negative of exponent used for astrometric RMS in linkage quality metric"
    )
    maxrms = lsst.pex.config.Field(
        dtype=float,
        default=1.0e5,
        doc="Maximum scaled RMS in km for a viable linkage"
    )
    max_oop = lsst.pex.config.Field(
        dtype=float,
        default=10000.0,
        doc="Maximum scaled out-of-plane RMS in km for a viable linkage"
    )
    rejfrac = lsst.pex.config.Field(
        dtype=float,
        default=0.5,
        doc="Maximum fraction of points that can be rejected"
    )
    maxrejnum = lsst.pex.config.Field(
        dtype=int,
        default=50,
        doc="Maximum number of points that can be rejected"
    )
    max_astrom_rms = lsst.pex.config.Field(
        dtype=float,
        default=1.0,
        doc="Maximum astrometric RMS in arcsec relative to best-fit orbit"
    )
    minobsnights = lsst.pex.config.Field(
        dtype=int,
        default=3,
        doc="Minimum number of distinct observing nights for a valid linkage"
    )
    minpointnum = lsst.pex.config.Field(
        dtype=int,
        default=6,
        doc="Minimum number of individual detections for a valid linkage"
    )
    use_heliovane = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Are we analyzing data from heliovane rather than heliolinc?"
    )
    verbose = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Prints monitoring output."
    )
    doLinkPlanarity = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Whether to use new linkPlanarity method, which approaches full completeness with faster runtime"
    )
    ecc_penalty = lsst.pex.config.Field(
        dtype=float,
        default=1.5,
        doc="Penalty factor for high-eccentricity orbits in linkage quality metric"
    )


class LinkPurifyTask(lsst.pipe.base.PipelineTask):
    ConfigClass = LinkPurifyConfig
    _DefaultName = "linkPurify"

    def run(self, sspVisitInputs, sspTrackletSources, sspBalancedLinkages, sspBalancedLinkageSources):
        """doc string 
           here
        """

        # copy all config parameters from the Task's config object
        # to heliolinx's native config object.

        if len(sspBalancedLinkages) == 0:
            raise NoWorkFound

        config = hl.LinkPurifyConfig()
        allvars = [item for item in dir(hl.LinkPurifyConfig) if not item.startswith("_")]
        for var in allvars:
            setattr(config, var, getattr(self.config, var))
        print('ecc', self.config.ecc_penalty, config.ecc_penalty)

        if self.config.doLinkPlanarity:
            (
                sspPurifiedLinkages, sspPurifiedLinkageSources
            ) = hl.linkPlanarity(config,
                                 utils.df2numpy(sspVisitInputs,      "hlimage"),
                                 utils.df2numpy(sspTrackletSources,  "hldet"),
                                 utils.df2numpy(sspBalancedLinkages,         "hlclust"),
                                 utils.df2numpy(sspBalancedLinkageSources,   "longpair"),
                                )

        else:
            (
                sspPurifiedLinkages, sspPurifiedLinkageSources
            ) = hl.linkPurify(config,
                              utils.df2numpy(sspVisitInputs,      "hlimage"),
                              utils.df2numpy(sspTrackletSources,  "hldet"),
                              utils.df2numpy(sspBalancedLinkages,         "hlclust"),
                              utils.df2numpy(sspBalancedLinkageSources,   "longpair"),
                             )

        return lsst.pipe.base.Struct(sspPurifiedLinkages=sspPurifiedLinkages,
                                     sspPurifiedLinkageSources=sspPurifiedLinkageSources
                                     )
