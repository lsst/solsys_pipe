#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ("LoadBalanceConfig", "LoadBalanceTask")

import logging

import astropy.table as tb
from lsst.pex.config import Field
from lsst.pipe.base import (
    NoWorkFound,
    PipelineTask,
    PipelineTaskConfig,
    PipelineTaskConnections,
    Struct,
    connectionTypes,
)

_LOG = logging.getLogger(__name__)


class LoadBalanceConnections(
    PipelineTaskConnections, dimensions=["instrument", "day_obs", "ssp_hypothesis_table"]
):
    sspLinkageList = connectionTypes.Input(
        doc="List of linkage tables to be load-balanced.",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"],
        storageClass="ArrowAstropy",
        name="ssp_linkages",
        multiple=True,
        deferLoad=True,
    )
    sspLinkageCounts = connectionTypes.Input(
        doc="Row counts of linkage tables to be load-balanced.",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"],
        storageClass="ArrowAstropy",
        name="ssp_linkages.rowcount",
        multiple=True,
    )
    sspLinkageSourceList = connectionTypes.Input(
        doc="List of linkage source tables corresponding to the linkage tables.",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"],
        storageClass="ArrowAstropy",
        name="ssp_linkage_sources",
        multiple=True,
        deferLoad=True,
    )
    sspLoadBalancedLinkages = connectionTypes.Output(
        doc="List of load-balanced linkage tables.",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_balanced_index"],
        storageClass="ArrowAstropy",
        name="ssp_balanced_linkages",
        multiple=True,
    )
    sspLoadBalancedLinkageSources = connectionTypes.Output(
        doc="List of load-balanced linkage source tables.",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_balanced_index"],
        storageClass="ArrowAstropy",
        name="ssp_balanced_linkage_sources",
        multiple=True,
    )


class LoadBalanceConfig(PipelineTaskConfig, pipelineConnections=LoadBalanceConnections):
    num_linkrefine_indices = Field(dtype=int, default=10, doc="Number of linkRefine quanta / output datasets")


class LoadBalanceTask(PipelineTask):
    ConfigClass = LoadBalanceConfig
    _DefaultName = "loadBalance"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        n_tables = len(inputs.sspLinkageList)
        if n_tables == 0:
            raise NoWorkFound

        _LOG.info(f"Concatenating {n_tables} linkage tables")

        n_ind = self.config.num_linkrefine_indices
        n_linkages = sum(inputs.sspLinkageCounts)
        n_target = int(n_linkages/n_ind) + 1
        leftoverLinkages, leftoverSources = None, None
        n_linkages_processed, i_output = 0, 0

        for i in range(n_tables):
            _LOG.info(f"Reading input table {i}")
            sspLinkage = inputs.sspLinkageList[i].get()
            sspLinkage["clusternum"] += n_linkages_processed
            sspLinkageSource = inputs.sspLinkageSourceList[i].get()
            sspLinkageSource["i1"] += n_linkages_processed
            n_linkages_processed += len(sspLinkage)
            if leftoverLinkages is None:
                leftoverLinkages = sspLinkage
                leftoverSources = sspLinkageSource
            else:
                leftoverLinkages = tb.vstack([leftoverLinkages, sspLinkage])
                leftoverSources = tb.vstack([leftoverSources, sspLinkageSource])

            while len(leftoverLinkages) >= n_target:

                outputLinkages = leftoverLinkages[:n_target]
                cutoff_clusternum = outputLinkages['clusternum'][-1]
                outputSources = leftoverSources[leftoverSources['i1'] <= cutoff_clusternum]

                leftoverLinkages = leftoverLinkages[n_target:]
                leftoverSources = leftoverSources[leftoverSources['i1'] > cutoff_clusternum]

                _LOG.info(f"Writing output tables {i_output}")
                linkageRef = outputRefs.sspLoadBalancedLinkages[i_output]
                butlerQC.put(outputLinkages, linkageRef)
                sourceRef = outputRefs.sspLoadBalancedLinkageSources[i_output]
                butlerQC.put(outputSources, sourceRef)

                i_output += 1
        if i_output < n_ind - 1:
            linkageRef = outputRefs.sspLoadBalancedLinkages[i_output]
            butlerQC.put(outputLinkages, linkageRef)
            sourceRef = outputRefs.sspLoadBalancedLinkageSources[i_output]
            butlerQC.put(outputSources, sourceRef)