import heliolinx.heliolinx as hl
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
import pandas as pd
import numpy as np
from . import utils

class MakeTrackletsConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=["instrument"]):
    sspDiaSourceInputs = connectionTypes.Input(
        doc="Table of unattributed sources",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name="sspDiaSourceInputs"
    )
    sspVisitInputs = connectionTypes.PrerequisiteInput(
        doc="visit stats plus observer coordinates",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name="sspVisitInputs"
    )
    sspTrackletSources = connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name="sspTrackletSources"
    )
    sspTracklets = connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name="sspTracklets"
    )
    sspTrackletToSource = connectionTypes.Output(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument"],
        storageClass="DataFrame",
        name="sspTrackletToSource"
    )


def getImageTimeTol():
    return 2./(24.*3600.)


class MakeTrackletsConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=MakeTrackletsConnections):
    mintrkpts = lsst.pex.config.Field(
        dtype=int,
        default=2,
        doc="minimum number of sources to qualify as a tracklet"
    )
    imagetimetol = lsst.pex.config.Field(
        dtype=float,
        default=getImageTimeTol(),
        doc="Tolerance for matching image time, in days: e.g. 1 second"
    )
    maxvel = lsst.pex.config.Field(
        dtype=float,
        default=1.5,
        doc="Default max angular velocity in deg/day."
    )
    minvel = lsst.pex.config.Field(
        dtype=float,
        default=0,
        doc="Min angular velocity in deg/day"
    )
    exptime = lsst.pex.config.Field(
        dtype=float,
        default=30,
        doc="FIXME WITH GOOD DOCUMENTATION. Exposure time"
    )
    minarc = lsst.pex.config.Field(
        dtype=float,
        default=0,
        doc="Min total angular arc in arcseconds."
    )
    maxtime = lsst.pex.config.Field(
        dtype=float,
        default=1.5/24,
        doc="Max inter-image time interval, in days."
    )
    mintime = lsst.pex.config.Field(
        dtype=float,
        default=1/86400,
        doc="Minimum inter-image time interval, in days."
    )
    imagerad = lsst.pex.config.Field(
        dtype=float,
        default=2.0,
        doc="radius from image center to most distant corner (deg)"
    )
    maxgcr = lsst.pex.config.Field(
        dtype=float,
        default=0.5,
        doc="Default maximum Great Circle Residual allowed for a valid tracklet (arcsec)"
    )
    timespan = lsst.pex.config.Field(
        dtype=float,
        default=14.0,
        doc="Default time to look back before most recent data (days)"
    )
    siglenscale = lsst.pex.config.Field(
        dtype=float,
        default=0.5,
        doc="????"
    )
    sigpascale = lsst.pex.config.Field(
        dtype=float,
        default=1.0,
        doc="????"
    )
    max_netl = lsst.pex.config.Field(
        dtype=int,
        default=2,
        doc="Maximum non-exclusive tracklet length"
    )
    forcerun = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Pushes through all but the immediately fatal errors."
    )
    verbose = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Prints monitoring output."
    )

class MakeTrackletsTask(lsst.pipe.base.PipelineTask):
    ConfigClass = MakeTrackletsConfig
    _DefaultName = "makeTracklets"

    def run(self, sspDiaSourceInputs, sspVisitInputs):
        """doc string
           here
        """

        # copy all config parameters from the Task's config object
        # to heliolinc's native config object.
        config = hl.MakeTrackletsConfig()
        allvars = [item for item in dir(hl.MakeTrackletsConfig) if not item.startswith("__")]
        for var in allvars:
            setattr(config, var, getattr(self.config, var))

        # convert dataframes to numpy array with dtypes that heliolinc expects
        (dets, tracklets, trac2det) = hl.makeTracklets(config,
                                                       utils.df2numpy(sspDiaSourceInputs, "hldet"),
                                                       utils.df2numpy(sspVisitInputs,     "hlimage"),
                                                      )

        # Do something about trailed sources
        return lsst.pipe.base.Struct(sspTrackletSources=dets,
                                     sspTracklets=tracklets,
                                     sspTrackletToSource=trac2det
                                     )
