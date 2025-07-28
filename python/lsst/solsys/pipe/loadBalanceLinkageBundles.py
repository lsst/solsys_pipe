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

import astropy.table
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
    )
    sspLinkageSourceList = connectionTypes.Input(
        doc="List of linkage source tables corresponding to the linkage tables.",
        dimensions=["day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"],
        storageClass="ArrowAstropy",
        name="ssp_linkage_sources",
        multiple=True,
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
        outputs = self.run(**inputs)
        n = len(outputs.sspLoadBalancedLinkageList)

        for i in range(n):
            dataId = outputRefs.sspLoadBalancedLinkages[i]
            butlerQC.put(outputs.sspLoadBalancedLinkageList[i], dataId)
            dataId = outputRefs.sspLoadBalancedLinkageSources[i]
            butlerQC.put(outputs.sspLoadBalancedLinkageSourceList[i], dataId)

    def run(self, sspLinkageList, sspLinkageSourceList):
        """Load balance the linkage tables across the specified number of
        indices.

        Parameters
        ----------
        sspLinkageList : `list` of `astropy.table.Table`
            List of linkage tables to be load-balanced.
        sspLinkageSourceList : `list` of `astropy.table.Table`
            List of linkage source tables corresponding to the linkage tables.

        Returns
        -------
        `Struct`
            A structure containing two lists:
            - ``sspLoadBalancedLinkageList``: List of load-balanced linkage
              tables.
            - ``sspLoadBalancedLinkageSourceList``: List of load-balanced
              linkage source tables.
        """
        if len(sspLinkageList) == 0:
            raise NoWorkFound

        n_ind = self.config.num_linkrefine_indices
        _LOG.info(f"Concatenating {len(sspLinkageList)} linkage tables")
        n = 0
        for i in range(len(sspLinkageList)):
            sspLinkage = sspLinkageList[i]
            sspLinkage["clusternum"] += n
            sspLinkageList[i] = sspLinkage
            sspLinkageSource = sspLinkageSourceList[i]
            sspLinkageSource["i1"] += n
            sspLinkageSourceList[i] = sspLinkageSource
            n += len(sspLinkage)
        sspLinkage = astropy.table.vstack(sspLinkageList)
        sspLinkageSource = astropy.table.vstack(sspLinkageSourceList)
        sspLinkage["loadBalanceIndex"] = sspLinkage["clusternum"] % n_ind
        sspLinkageSource["loadBalanceIndex"] = sspLinkageSource["i1"] % n_ind
        sspLoadBalancedLinkageList = [t for t in sspLinkage.group_by("loadBalanceIndex").groups]
        sspLoadBalancedLinkageSourceList = [t for t in sspLinkageSource.group_by("loadBalanceIndex").groups]
        for i in range(n_ind):
            if len(sspLoadBalancedLinkageList) <= i:
                sspLoadBalancedLinkageList.append(astropy.table.Table())
                sspLoadBalancedLinkageSourceList.append(astropy.table.Table())
            else:
                sspLoadBalancedLinkageList[i]["clusternum"] //= n_ind
                sspLoadBalancedLinkageSourceList[i]["i1"] //= n_ind
                sspLoadBalancedLinkageList[i].remove_column("loadBalanceIndex")
                sspLoadBalancedLinkageSourceList[i].remove_column("loadBalanceIndex")

        return Struct(
            sspLoadBalancedLinkageList=sspLoadBalancedLinkageList,
            sspLoadBalancedLinkageSourceList=sspLoadBalancedLinkageSourceList,
        )
