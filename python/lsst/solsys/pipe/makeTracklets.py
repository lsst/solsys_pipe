import heliolinx.heliolinx as hl
import lsst.pex.config
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
import pandas as pd
import numpy as np
from . import utils
import astropy.units as u
from datetime import datetime

import lsst.pex.config as pexConfig
from lsst.pipe.base import Struct, NoWorkFound
import lsst.pipe.base.connectionTypes as connTypes
from lsst.verify import Measurement, Datum
from lsst.verify.tasks import AbstractMetadataMetricTask, MetricTask, MetricComputationError

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pipe.base.connectionTypes as cT
import numpy as np
from lsst.pex.config import Config, ConfigField, ConfigurableField, Field
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.tasks.background import (
    FocalPlaneBackground,
    FocalPlaneBackgroundConfig,
    MaskObjectsTask,
    SkyMeasurementTask,
)

from collections import defaultdict
import dataclasses
import functools
import logging
import numbers
import os

import numpy as np
import pandas as pd
import astropy.table
from astro_metadata_translator.headers import merge_headers

import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
from lsst.daf.butler.formatters.parquet import pandas_to_astropy
from lsst.pipe.base import NoWorkFound, connectionTypes
import lsst.afw.table as afwTable
from lsst.afw.image import ExposureSummaryStats, ExposureF
from lsst.meas.base import SingleFrameMeasurementTask, DetectorVisitIdGeneratorConfig
from lsst.obs.base.utils import strip_provenance_from_fits_header

from lsst.pipe.tasks.postprocess import TableVStack


class MakeTrackletsConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=["instrument", "day_obs"]):
    sspDiaSourceInputs = connectionTypes.Input(
        doc="Table of unattributed sources",
        dimensions=["instrument", "day_obs"],
        storageClass="DataFrame",
        name="dia_source_dayobs"
    )
    sspVisitInputs = connectionTypes.PrerequisiteInput(
        doc="visit stats plus observer coordinates",
        dimensions=["instrument", "day_obs"],
        storageClass="DataFrame",
        name="sspVisitInputs"
    )
    sspTrackletSources = connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="sspTrackletSources"
    )
    sspTracklets = connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="sspTracklets"
    )
    sspTrackletToSource = connectionTypes.Output(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="sspTrackletToSource"
    )
    sspVisitHeliolincInputs = connectionTypes.Output(
        doc="visit stats plus observer coordinates formatted for heliolinc",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="sspVisitHeliolincInputs",
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
        doc="Default exposure time (overriden by sspVisitInputs)"
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
    siglenscale = lsst.pex.config.Field(
        dtype=float,
        default=0.5,
        doc="Default scaling from trail length to trail length uncertainty"
    )
    sigpascale = lsst.pex.config.Field(
        dtype=float,
        default=1.0,
        doc="Default scaling from trail length to trail angle uncertainty"
    )
    max_netl = lsst.pex.config.Field(
        dtype=int,
        default=2,
        doc="Maximum non-exclusive tracklet length"
    )
    verbose = lsst.pex.config.Field(
        dtype=int,
        default=0,
        doc="Prints monitoring output."
    )
    time_offset = lsst.pex.config.Field(
        dtype=float,
        default=0,
        doc="Offset in seconds to change timescale (TAI to UTC, for example)"
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
