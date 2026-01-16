import heliolinx.heliolinx as hl
import heliolinx.solarsyst_dyn_geo as sdg
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
from . import utils
import astropy.table as tb

# This task purifies and de-duplicates candidate asteroid linkages output
# by heliolinc

class LinkMergeConnections(lsst.pipe.base.PipelineTaskConnections,
                             dimensions=("instrument", "day_obs",
                                         "ssp_hypothesis_table")):
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
    sspPurifiedLinkages = connectionTypes.Input(
        doc="one line summary of each linkage",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_balanced_index", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_purified_linkages",
        multiple=True
    )
    sspPurifiedLinkageSources = connectionTypes.Input(
        doc="indices connecting linkages to trackletSources",
        dimensions=("instrument", "ssp_hypothesis_table", "ssp_balanced_index", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_purified_linkage_sources",
        multiple=True
    )
    sspMergedLinkages = connectionTypes.Output(
        doc="one line summary of each linkage",
        dimensions=("instrument", "ssp_hypothesis_table", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_merged_linkages",
    )
    sspMergedLinkageSources = connectionTypes.Output(
        doc="indices connecting linkages to trackletSources",
        dimensions=("instrument", "ssp_hypothesis_table", "day_obs"),
        storageClass="ArrowAstropy",
        name = "ssp_merged_linkage_sources",
    )


class LinkMergeConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=LinkMergeConnections):
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
        doc="Exponent used for number of points in linkage quality metric"
    )
    nightpow = lsst.pex.config.Field(
        dtype=int,
        default=-1,
        doc="Exponent used for number of nights in linkage quality metric"
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
        default=1000.0,
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
    ecc_penalty = lsst.pex.config.Field(
        dtype=float,
        default=1.5,
        doc="Penalty factor for high-eccentricity orbits in linkage quality metric"
    ) 


class LinkMergeTask(lsst.pipe.base.PipelineTask):
    ConfigClass = LinkMergeConfig
    _DefaultName = "linkMerge"

    def run(self, sspVisitInputs, sspTrackletSources, sspPurifiedLinkages, sspPurifiedLinkageSources):
        """doc string 
           here
        """

        # consolidate stuff
        n = 0
        n_tables = len(sspPurifiedLinkages)
        for i in range(n_tables):
            sspLinkage = sspPurifiedLinkages[i]
            sspLinkage['clusternum'] += n
            sspPurifiedLinkages[i] = sspLinkage
            sspLinkageSource = sspPurifiedLinkageSources[i]
            sspLinkageSource['i1'] += n
            sspPurifiedLinkageSources[i] = sspLinkageSource
            n += len(sspLinkage)

        sspPurifiedLinkages = tb.vstack(sspPurifiedLinkages)
        sspPurifiedLinkageSources = tb.vstack(sspPurifiedLinkageSources)

        # copy all config parameters from the Task's config object
        # to heliolinx's native config object.
        config = hl.LinkPurifyConfig()
        allvars = [item for item in dir(hl.LinkPurifyConfig) if not item.startswith("_")]
        for var in allvars:
            setattr(config, var, getattr(self.config, var))

        (
            sspMergedLinkages, sspMergedLinkageSourceIndices
        ) = hl.linkPurify(config,
                          utils.df2numpy(sspVisitInputs,      "hlimage"),
                          utils.df2numpy(sspTrackletSources,  "hldet"),
                          utils.df2numpy(sspPurifiedLinkages,         "hlclust"),
                          utils.df2numpy(sspPurifiedLinkageSources,   "longpair"),
                         )

        sspMergedLinkageSources = sspTrackletSources[sspMergedLinkageSourceIndices['i2']]
        sspMergedLinkageSources['clusternum'] = sspMergedLinkageSourceIndices['i1']

        return lsst.pipe.base.Struct(sspMergedLinkages=sspMergedLinkages,
                                     sspMergedLinkageSources=sspMergedLinkageSources,
                                     )
