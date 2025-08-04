import heliolinx.heliolinx as hl
import lsst.pipe.base
from lsst.pipe.base import connectionTypes
import pandas as pd
import numpy as np
from . import utils
import astropy.units as u
from datetime import datetime

import lsst.pipe.base.connectionTypes as connTypes
from lsst.verify import Measurement, Datum
from lsst.verify.tasks import AbstractMetadataMetricTask, MetricTask, MetricComputationError

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pipe.base.connectionTypes as cT
from lsst.pex.config import Config, ConfigField, ConfigurableField, Field
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, NoWorkFound
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
import astropy.table
from astro_metadata_translator.headers import merge_headers

import lsst.geom
import lsst.pipe.base as pipeBase
import lsst.daf.base as dafBase
from lsst.daf.butler.formatters.parquet import pandas_to_astropy
import lsst.afw.table as afwTable
from lsst.afw.image import ExposureSummaryStats, ExposureF
from lsst.meas.base import SingleFrameMeasurementTask, DetectorVisitIdGeneratorConfig
from lsst.obs.base.utils import strip_provenance_from_fits_header

from lsst.pipe.tasks.postprocess import TableVStack

_LOG = logging.getLogger(__name__)
visitSummaryColumnRenameDict = {'MJD': 'MJD', 'boresightRa': 'RA', 'boresightDec': 'Dec', 'exposureTime': 'exptime'}

class MakeTrackletsConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=["instrument", "day_obs", "ssp_hypothesis_table"]):
    sspDiaSourceInputs = connectionTypes.Input(
        doc="Table of unattributed sources",
        dimensions=["instrument", "day_obs"],
        storageClass="DataFrame",
        name="dia_source_dayobs_14"
    )
    sspVisitInputs = connectionTypes.Input(
        doc="visit stats plus observer coordinates",
        dimensions=["instrument", "day_obs"],
        storageClass="DataFrame",
        name="visit_summary_dayobs_14"
    )
    sspTrackletSources = connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_sources"
    )
    sspTracklets = connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklets"
    )
    sspTrackletToSource = connectionTypes.Output(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_to_source"
    )
    sspVisitHeliolincInputs = connectionTypes.Output(
        doc="visit stats plus observer coordinates formatted for heliolinc",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_visit_heliolinc_inputs",
    )


diaSourceColumnRenameDict = {'diaSourceId': 'idstring', 'visit': 'image', 'midpointMjdTai': 'MJD',
                             'ra': 'RA', 'dec': 'Dec', 'trailLength': 'trail_len', 'trailAngle': 'trail_PA'}
visitSummaryColumnRenameDict = {'MJD': 'MJD', 'boresightRa': 'RA', 'boresightDec': 'Dec', 'exposureTime': 'exptime'}

class MakeTrackletsConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=MakeTrackletsConnections):
    mintrkpts = Field(
        dtype=int,
        default=2,
        doc="minimum number of sources to qualify as a tracklet"
    )
    imagetimetol = Field(
        dtype=float,
        default=1./(24.*3600.),
        doc="Tolerance for matching image time, in days (default: 1 second)"
    )
    maxvel = Field(
        dtype=float,
        default=1.5,
        doc="Default max angular velocity in deg/day."
    )
    minvel = Field(
        dtype=float,
        default=0,
        doc="Min angular velocity in deg/day"
    )
    exptime = Field(
        dtype=float,
        default=30,
        doc="Default exposure time in seconds (overriden by sspVisitInputs)"
    )
    minarc = Field(
        dtype=float,
        default=0,
        doc="Min total angular arc in arcseconds."
    )
    maxtime = Field(
        dtype=float,
        default=1.5/24,
        doc="Max inter-image time interval, in days."
    )
    mintime = Field(
        dtype=float,
        default=1/86400,
        doc="Minimum inter-image time interval, in days."
    )
    imagerad = Field(
        dtype=float,
        default=2.0,
        doc="radius from image center to most distant corner (deg)"
    )
    maxgcr = Field(
        dtype=float,
        default=0.5,
        doc="Default maximum Great Circle Residual allowed for a valid tracklet (arcsec)"
    )
    siglenscale = Field(
        dtype=float,
        default=0.5,
        doc="Default scaling from trail length to trail length uncertainty"
    )
    sigpascale = Field(
        dtype=float,
        default=1.0,
        doc="Default scaling from trail length to trail angle uncertainty"
    )
    max_netl = Field(
        dtype=int,
        default=2,
        doc="Maximum non-exclusive tracklet length"
    )
    verbose = Field(
        dtype=int,
        default=0,
        doc="Prints monitoring output."
    )
    time_offset = Field(
        dtype=float,
        default=0,
        doc="Offset in seconds to change timescale (TAI to UTC, for example)"
    )
    trkfrac = Field(
        dtype=float,
        default=0.3,
        doc="Minimum tracklet length in high-depth fields as a factor of depth."
    )
    matchrad = Field(
        dtype=float,
        default=0.2,
        doc="Radius of field matching in degrees"
    )
    use_lowmem = Field(
        dtype=bool,
        default=True,
        doc="Use the new low-memory makeTracklets algorithm. Config left to allow backwards-compatibility just in case."
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
        configs_to_transfer = ['mintrkpts', 'imagetimetol', 'maxvel', 'minvel', 'exptime', 'minarc',
                               'maxtime', 'mintime', 'imagerad', 'matchrad', 'trkfrac', 'maxgcr', 'siglenscale', 'sigpascale',
                               'max_netl', 'use_lowmem', 'verbose']  # should this have time_offset?
        for config_name in configs_to_transfer:
            setattr(config, config_name, getattr(self.config, config_name))

        hldet = utils.make_hldet(sspDiaSourceInputs)

        print(sspVisitInputs)
        (dets, tracklets, trac2det) = hl.makeTracklets(config,
                                                       hldet,
                                                       utils.make_hlimage(sspVisitInputs),
                                                      )
        _LOG.info(f'makeTracklets finished: {len(dets)} detections, {len(tracklets)} tracklets')
        # Do something about trailed sources
        return lsst.pipe.base.Struct(sspTrackletSources=dets,
                                     sspTracklets=tracklets,
                                     sspTrackletToSource=trac2det
                                     )
