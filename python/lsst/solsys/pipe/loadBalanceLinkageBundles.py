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

class LoadBalanceConnections(lsst.pipe.base.PipelineTaskConnections,
                               dimensions=["instrument", "day_obs", "ssp_hypothesis_table"]):
    sspLinkageList = connectionTypes.Input(
        doc="",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"],
        storageClass="ArrowAstropy",
        name="ssp_linkages",
        multiple=True,
    )
    sspLinkageSourceList = connectionTypes.Input(
        doc="",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"],
        storageClass="ArrowAstropy",
        name="ssp_linkage_sources",
        multiple=True,
    )
    sspLoadBalancedLinkages = connectionTypes.Output(
        doc="",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_balanced_index"],
        storageClass="ArrowAstropy",
        name="ssp_balanced_linkages",
        multiple=True,
    )
    sspLoadBalancedLinkageSources = connectionTypes.Output(
        doc="",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_balanced_index"],
        storageClass="ArrowAstropy",
        name="ssp_balanced_linkage_sources",
        multiple=True,
    )


class LoadBalanceConfig(lsst.pipe.base.PipelineTaskConfig, pipelineConnections=LoadBalanceConnections):
    num_linkrefine_indices = Field(
        dtype=int,
        default=10,
        doc="Number of linkRefine quanta / output datasets"
    )

class LoadBalanceTask(lsst.pipe.base.PipelineTask):
    ConfigClass = LoadBalanceConfig
    _DefaultName = "loadBalance"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(**inputs)
        n = len(outputs.sspLoadBalancedLinkageList)
        for i in range(n):
            dataId = outputRefs.sspLoadBalancedLinkages[i]
            butlerQC.put(outputs.sspLoadBalancedLinkageList[i], dataId)
            dataId = outputRefs.sspLoadBalancedLinkageSources[i]
            butlerQC.put(outputs.sspLoadBalancedLinkageSourceList[i], dataId)

    def run(self, sspLinkageList, sspLinkageSourceList):
        """doc string
           here
        """
        if len(sspLinkageList) == 0:
            raise NoWorkFound

        n_ind = self.config.num_linkrefine_indices
        _LOG.info(f'Concatenating {len(sspLinkageList)} linkage tables')
        n = 0
        for i in range(len(sspLinkageList)):
            sspLinkage = sspLinkageList[i]
            sspLinkage['clusternum'] += n
            sspLinkageList[i] = sspLinkage
            sspLinkageSource = sspLinkageSourceList[i]
            sspLinkageSource['i1'] += n
            sspLinkageSourceList[i] = sspLinkageSource
            n += len(sspLinkage)
        sspLinkage = astropy.table.vstack(sspLinkageList)
        sspLinkageSource = astropy.table.vstack(sspLinkageSourceList)
        sspLinkage['loadBalanceIndex'] = sspLinkage['clusternum'] % n_ind
        sspLinkageSource['loadBalanceIndex'] = sspLinkageSource['i1'] % n_ind
        sspLoadBalancedLinkageList = [t for t in sspLinkage.group_by('loadBalanceIndex').groups]
        sspLoadBalancedLinkageSourceList = [t for t in sspLinkageSource.group_by('loadBalanceIndex').groups]
        for i in range(n_ind):
            sspLoadBalancedLinkageList[i]['clusternum'] //= n_ind
            sspLoadBalancedLinkageSourceList[i]['i1'] //= n_ind
            sspLoadBalancedLinkageList[i].remove_column('loadBalanceIndex')
            sspLoadBalancedLinkageSourceList[i].remove_column('loadBalanceIndex')

        return lsst.pipe.base.Struct(sspLoadBalancedLinkageList = sspLoadBalancedLinkageList,
                                     sspLoadBalancedLinkageSourceList = sspLoadBalancedLinkageSourceList
                                     )
