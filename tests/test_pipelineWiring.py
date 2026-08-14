# This file is part of solsys_pipe.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Wiring checks for dataset names passed between SSP pipeline tasks."""

import unittest
from typing import NamedTuple

import lsst.utils.tests
from lsst.solsys.pipe.consolidateTracklets import (
    ConsolidateTrackletsConfig,
    ConsolidateTrackletsConnections,
)
from lsst.solsys.pipe.heliolinc import HeliolincConfig, HeliolincConnections
from lsst.solsys.pipe.linkMerge import LinkMergeConfig, LinkMergeConnections
from lsst.solsys.pipe.linkPurify import LinkPurifyConfig, LinkPurifyConnections
from lsst.solsys.pipe.loadBalanceLinkageBundles import (
    LoadBalanceConfig,
    LoadBalanceConnections,
)
from lsst.solsys.pipe.makeTracklets import (
    MakeTrackletsConfig,
    MakeTrackletsConnections,
)


class Wire(NamedTuple):
    """One table handed from an upstream task to a downstream task."""

    producer: str  # Task connection that writes a table
    consumer: str  # Task connection that reads that same table
    purpose: str  # Why is the table needed downstream?


# These are the pipeline's internal handoffs. If a task starts writing or
# reading a new internal table, add it here so the handoff is visible.
# Stage names match the task labels in the DIA pipeline YAML.
PIPELINE_WIRES = (
    Wire(
        "makeTracklets.outputVisitInfo",
        "consolidateTracklets.inputVisitSummaries",
        "Consolidate stacks the daily visit table MakeTracklets writes.",
    ),
    Wire(
        "makeTracklets.sspTrackletSources",
        "consolidateTracklets.inputTrackletSources",
        "Consolidate stacks the daily tracklet-source table MakeTracklets writes.",
    ),
    Wire(
        "makeTracklets.sspTracklets",
        "consolidateTracklets.inputTracklets",
        "Consolidate stacks the daily tracklet table MakeTracklets writes.",
    ),
    Wire(
        "makeTracklets.sspTrackletToSource",
        "consolidateTracklets.inputTrackletToSource",
        "Consolidate stacks the daily tracklet-to-source map MakeTracklets writes.",
    ),
    Wire(
        "consolidateTracklets.outputVisitSummaries",
        "heliolinc.sspVisitInputs",
        "Heliolinc needs the same 14-day visit table Consolidate writes.",
    ),
    Wire(
        "consolidateTracklets.outputTrackletSources",
        "heliolinc.sspTrackletSources",
        "Heliolinc needs the same 14-day tracklet-source table Consolidate writes.",
    ),
    Wire(
        "consolidateTracklets.outputTracklets",
        "heliolinc.sspTracklets",
        "Heliolinc needs the same 14-day tracklet table Consolidate writes.",
    ),
    Wire(
        "consolidateTracklets.outputTrackletToSource",
        "heliolinc.sspTrackletToSource",
        "Heliolinc needs the same tracklet-to-source map Consolidate writes.",
    ),
    Wire(
        "heliolinc.sspLinkage",
        "loadBalance.sspLinkageList",
        "LoadBalance splits the linkage tables Heliolinc writes.",
    ),
    Wire(
        "heliolinc.sspLinkageSources",
        "loadBalance.sspLinkageSourceList",
        "LoadBalance keeps linkage-source rows paired with those linkages.",
    ),
    Wire(
        "loadBalance.sspLoadBalancedLinkages",
        "linkPurify.sspBalancedLinkages",
        "LinkPurify cleans the balanced linkage partitions.",
    ),
    Wire(
        "loadBalance.sspLoadBalancedLinkageSources",
        "linkPurify.sspBalancedLinkageSources",
        "LinkPurify needs the source rows for each balanced linkage partition.",
    ),
    Wire(
        "linkPurify.sspPurifiedLinkages",
        "linkMerge.sspPurifiedLinkages",
        "LinkMerge combines all purified linkage partitions.",
    ),
    Wire(
        "linkPurify.sspPurifiedLinkageSources",
        "linkMerge.sspPurifiedLinkageSources",
        "LinkMerge combines the source rows for those purified linkages.",
    ),
    Wire(
        "consolidateTracklets.outputVisitSummaries",
        "linkPurify.sspVisitInputs",
        "LinkPurify also needs the 14-day visit table for orbit checks.",
    ),
    Wire(
        "consolidateTracklets.outputTrackletSources",
        "linkPurify.sspTrackletSources",
        "LinkPurify also needs the 14-day tracklet-source table.",
    ),
    Wire(
        "consolidateTracklets.outputVisitSummaries",
        "linkMerge.sspVisitInputs",
        "LinkMerge also needs the 14-day visit table for orbit checks.",
    ),
    Wire(
        "consolidateTracklets.outputTrackletSources",
        "linkMerge.sspTrackletSources",
        "LinkMerge also needs the 14-day tracklet-source table.",
    ),
)

# LinkMerge is the final stage in this test, so its outputs do not have to
# feed another in-repo task here. Every earlier output should be consumed
# unless listed as a side output below.
NON_FINAL_STAGES = (
    "makeTracklets",
    "consolidateTracklets",
    "heliolinc",
    "loadBalance",
    "linkPurify",
)

# MakeTracklets is the first stage in this test. Its inputs come from the DIA
# source side, so they are allowed to be outside this map.
FIRST_STAGE = "makeTracklets"

# Daily DIA sources are still a useful side product, but no later task in this
# SSP linkage chain reads them.
SIDE_OUTPUTS = {
    "makeTracklets": {"outputDiaTable"},
}

# Row count is a component of Heliolinc's linkage table, not a separate task
# output. It still counts as explained pipeline wiring.
SPECIAL_INTERNAL_INPUTS = {
    "loadBalance": {"sspLinkageCounts"},
}


class PipelineWiringTestCase(lsst.utils.tests.TestCase):
    """Check that pipeline task outputs feed the right downstream inputs."""

    def setUp(self):
        make_tracklets_config = MakeTrackletsConfig()
        make_tracklets_config.consolidateVisitTables = False
        self.connection_sets = {
            "makeTracklets": MakeTrackletsConnections(config=make_tracklets_config),
            "consolidateTracklets": ConsolidateTrackletsConnections(
                config=ConsolidateTrackletsConfig()
            ),
            "heliolinc": HeliolincConnections(config=HeliolincConfig()),
            "loadBalance": LoadBalanceConnections(config=LoadBalanceConfig()),
            "linkPurify": LinkPurifyConnections(config=LinkPurifyConfig()),
            "linkMerge": LinkMergeConnections(config=LinkMergeConfig()),
        }

    def _split_endpoint(self, endpoint):
        return endpoint.split(".", 1)

    def _connection(self, endpoint):
        stage, connection_name = self._split_endpoint(endpoint)
        return self.connection_sets[stage].allConnections[connection_name]

    def _connection_names_used_as(self, endpoint_attr):
        used = {}
        for wire in PIPELINE_WIRES:
            stage, connection_name = self._split_endpoint(
                getattr(wire, endpoint_attr)
            )
            used.setdefault(stage, set()).add(connection_name)
        return used

    def test_each_wire_reads_the_table_written_upstream(self):
        for wire in PIPELINE_WIRES:
            output = self._connection(wire.producer)
            input_ = self._connection(wire.consumer)
            label = f"{wire.producer} -> {wire.consumer}"
            self.assertEqual(
                output.name,
                input_.name,
                msg=f"{label}: {wire.purpose}",
            )
            self.assertEqual(
                output.storageClass,
                input_.storageClass,
                msg=f"{label}: storage class must match for the same table.",
            )

    def test_pipeline_handoffs_are_astropy_tables(self):
        for wire in PIPELINE_WIRES:
            output = self._connection(wire.producer)
            input_ = self._connection(wire.consumer)
            label = f"{wire.producer} -> {wire.consumer}"
            self.assertEqual(output.storageClass, "ArrowAstropy", msg=label)
            self.assertEqual(input_.storageClass, "ArrowAstropy", msg=label)

    def test_no_upstream_output_is_left_unread(self):
        used_outputs = self._connection_names_used_as("producer")
        for stage in NON_FINAL_STAGES:
            expected = set(used_outputs.get(stage, set()))
            expected.update(SIDE_OUTPUTS.get(stage, set()))
            self.assertEqual(
                set(self.connection_sets[stage].outputs),
                expected,
                msg=f"{stage} has an output not listed in PIPELINE_WIRES.",
            )

    def test_no_internal_input_is_unexplained(self):
        used_inputs = self._connection_names_used_as("consumer")
        for stage, special_inputs in SPECIAL_INTERNAL_INPUTS.items():
            used_inputs.setdefault(stage, set()).update(special_inputs)

        for stage, connections in self.connection_sets.items():
            if stage == FIRST_STAGE:
                continue
            self.assertEqual(
                set(connections.inputs),
                used_inputs.get(stage, set()),
                msg=f"{stage} has an input not listed in PIPELINE_WIRES.",
            )

    def test_make_tracklets_is_in_dia_mode(self):
        # The DIA pipeline sets consolidateVisitTables=False. In that mode
        # MakeTracklets starts from DIA sources and does not ask the pipeline
        # for a separate visit-summary input.
        self.assertNotIn(
            "inputVisitSummaries",
            self.connection_sets["makeTracklets"].inputs,
        )
        self.assertIn("inputDiaTables", self.connection_sets["makeTracklets"].inputs)

    def test_load_balance_rowcount_uses_heliolinc_linkage_component(self):
        self.assertEqual(
            self._connection("loadBalance.sspLinkageCounts").name,
            self._connection("heliolinc.sspLinkage").name + ".rowcount",
        )
        self.assertEqual(
            self._connection("loadBalance.sspLinkageCounts").storageClass,
            "int",
        )


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
