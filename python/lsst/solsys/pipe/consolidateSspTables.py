# This file is part of solsys_pipe
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""A per-dayobs task which takes `goodSeeingDiff_assocDiaSrc` or
`goodSeeingDiff_diaSrcTable` from the last 14 (config-settable) days and
consolidates them into a `sspVisitInputs` as defined in `makeTracklets.py`
"""

__all__ = [
    "ConsolidateSspTablesTask", "ConsolidateSspTablesConfig",
]

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

# from .functors import CompositeFunctor, Column

from lsst.pipe.tasks.postprocess import TableVStack

# From pipe_tasks/postprocessing!

"""A per-dayobs task which takes `goodSeeingDiff_assocDiaSrc` or
`goodSeeingDiff_diaSrcTable` from the last 14 (config-settable) days and
consolidates them into a `sspVisitInputs` as defined in `makeTracklets.py`
"""


class ConsolidateSspTablesConnections(PipelineTaskConnections,
                                            defaultTemplates={"coaddName": "deepSeeing"},
                                            dimensions=("instrument", "day_obs")):
    # connections.inputCatalogs: goodSeeingDiff_diaSrcTable
    # inputCatalogs = connectionTypes.Input(
    #     doc="Input per-detector Source Tables",
    #     name="{catalogType}sourceTable",
    #     storageClass="ArrowAstropy",
    #     dimensions=("instrument", "day_obs", "detector"),
    #     multiple=True,
    #     deferLoad=True,
    # )
    inputCatalogs = connTypes.Input(
        doc="Our input Source Tables to be concatenated",
        name="dia_source_detector",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        multiple=True,
    )
    outputCatalog = connectionTypes.Output(
        doc="Per-dayobs concatenation of Source Table",
        name="dia_source_dayobs",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "day_obs")
    )


class ConsolidateSspTablesConfig(pipeBase.PipelineTaskConfig,
                                       pipelineConnections=ConsolidateSspTablesConnections):
    pass


class ConsolidateSspTablesTask(pipeBase.PipelineTask):
    """Concatenate `sourceTable` list into a per-visit `sourceTable_visit`
    """
    _DefaultName = "consolidateDiaSourceTable"
    ConfigClass = ConsolidateSspTablesConfig

    inputDataset = "sourceTable"
    outputDataset = "sourceTable_visit"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        from .makeWarp import reorderRefs

        detectorOrder = [ref.dataId["detector"] for ref in inputRefs.inputCatalogs]
        detectorOrder.sort()
        inputRefs = reorderRefs(inputRefs, detectorOrder, dataIdKey="day_obs")
        inputs = butlerQC.get(inputRefs)
        self.log.info("Concatenating %s per-dayobs Source Tables",
                      len(inputs["inputCatalogs"]))
        table = TableVStack.vstack_handles(inputs["inputCatalogs"])
        butlerQC.put(pipeBase.Struct(outputCatalog=table), outputRefs)

