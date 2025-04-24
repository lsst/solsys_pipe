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
    sspVisitInputs = connectionTypes.Input(
        doc="visit stats plus observer coordinates",
        dimensions=["instrument", "day_obs"],
        storageClass="DataFrame",
        name="visit_summary_dayobs"
    )
    sspTrackletSources = connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_sources"
    )
    sspTracklets = connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_tracklets"
    )
    sspTrackletToSource = connectionTypes.Output(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_to_source"
    )
    sspVisitHeliolincInputs = connectionTypes.Output(
        doc="visit stats plus observer coordinates formatted for heliolinc",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_visit_heliolinc_inputs",
    )


def getImageTimeTol():
    return 2./(24.*3600.)

diaSourceColumnRenameDict = {'diaSourceId': 'idstring', 'visit': 'image', 'midpointMjdTai': 'MJD',
                             'ra': 'RA', 'dec': 'Dec', 'trailLength': 'trail_len', 'trailAngle': 'trail_PA'}
""" # Missing
float mag
float sigmag
char obscode[MINSTRINGLEN]
long known_obj
long det_qual"""
# ['diaSourceId', 'detector', 'band', 'diaObjectId', 'ssObjectId', 'parentDiaSourceId', 'midpointMjdTai', 'bboxSize', 'time_processed', 'ra', 'dec', 'raErr', 'decErr', 'ra_dec_Cov', 'x', 'y', 'xErr', 'yErr', 'apFlux', 'apFluxErr', 'snr', 'psfFlux', 'psfFluxErr', 'psfChi2', 'psfNdata', 'trailFlux', 'trailRa', 'trailDec', 'trailLength', 'trailAngle', 'dipoleMeanFlux', 'dipoleMeanFluxErr', 'dipoleFluxDiff', 'dipoleFluxDiffErr', 'dipoleLength', 'dipoleAngle', 'dipoleChi2', 'isDipole', 'dipoleFitAttempted', 'dipoleNdata', 'scienceFlux', 'scienceFluxErr', 'ixx', 'iyy', 'ixy', 'ixxPSF', 'iyyPSF', 'ixyPSF', 'extendedness', 'reliability', 'pixelFlags', 'pixelFlags_offimage', 'pixelFlags_edge', 'pixelFlags_interpolated', 'pixelFlags_saturated', 'pixelFlags_cr', 'pixelFlags_bad', 'pixelFlags_suspect', 'pixelFlags_interpolatedCenter', 'pixelFlags_saturatedCenter', 'pixelFlags_crCenter', 'pixelFlags_suspectCenter', 'centroid_flag', 'apFlux_flag', 'apFlux_flag_apertureTruncated', 'psfFlux_flag', 'psfFlux_flag_noGoodPixels', 'psfFlux_flag_edge', 'forced_PsfFlux_flag', 'forced_PsfFlux_flag_noGoodPixels', 'forced_PsfFlux_flag_edge', 'shape_flag', 'shape_flag_no_pixels', 'shape_flag_not_contained', 'shape_flag_parent_source', 'trail_flag_edge', 'pixelFlags_streak', 'pixelFlags_streakCenter', 'pixelFlags_injected', 'pixelFlags_injectedCenter', 'pixelFlags_injected_template', 'pixelFlags_injected_templateCenter', 'pixelFlags_nodata', 'pixelFlags_nodataCenter']
visitSummaryColumnRenameDict = {'MJD': 'MJD', 'boresightRa': 'RA', 'boresightDec': 'Dec', 'exposureTime': 'exptime'}
"""   double MJD;
  double RA;
  double Dec;
  char obscode[MINSTRINGLEN];
  double X;
  double Y;
  double Z;
  double VX;
  double VY;
  double VZ;
  long startind;
  long endind;
  double exptime; // Exposure time in seconds"""

class MakeTrackletsConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=MakeTrackletsConnections):
    obscode = lsst.pex.config.Field(
        dtype=str,
        default="X05",
        doc="MPC code of observatory"
    )
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
        configs_to_transfer = ['mintrkpts', 'imagetimetol', 'maxvel', 'minvel', 'exptime', 'minarc',
                               'maxtime', 'mintime', 'imagerad', 'maxgcr', 'siglenscale', 'sigpascale',
                               'max_netl', 'verbose'] # should this have time_offset?
        for config_name in configs_to_transfer:
            setattr(config, config_name, getattr(self.config, config_name))

        sspDiaSourceInputs = sspDiaSourceInputs.rename(columns = diaSourceColumnRenameDict)
        sspDiaSourceInputs = sspDiaSourceInputs[['MJD', 'RA', 'Dec', 'idstring']]
        sspDiaSourceInputs['idstring'] = sspDiaSourceInputs['idstring'].astype(str)
        sspDiaSourceInputs['obscode'] = self.config.obscode

        sspVisitInputs = sspVisitInputs.rename(columns = visitSummaryColumnRenameDict)
        sspVisitInputs = sspVisitInputs[['MJD', 'RA', 'Dec', 'exptime']]
        sspVisitInputs['obscode'] = self.config.obscode
        # convert dataframes to numpy array with dtypes that heliolinc expects
        (dets, tracklets, trac2det) = hl.makeTracklets(config,
                                                       utils.make_hldet(sspDiaSourceInputs),
                                                       utils.make_hlimage(sspVisitInputs),
                                                      )

        # Do something about trailed sources
        return lsst.pipe.base.Struct(sspTrackletSources=dets,
                                     sspTracklets=tracklets,
                                     sspTrackletToSource=trac2det
                                     )
