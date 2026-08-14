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

"""Config and connection checks for LinkMerge.

LinkMerge shares a lot of knobs with LinkPurify. If they drift apart,
the merged catalog gets filtered with different thresholds than the
per-bundle output it came from. The shared-knob test below pins that.
"""

import unittest

import lsst.pipe.base as pipeBase
import lsst.utils.tests
from lsst.solsys.pipe.linkMerge import (
    LinkMergeConfig,
    LinkMergeConnections,
    LinkMergeTask,
)
from lsst.solsys.pipe.linkPurify import LinkPurifyConfig

from synthetic import EXPECTED_LINKAGE_COLS, make_linkage_inputs


SHARED_KNOBS = (
    "minobsnights",
    "minpointnum",
    "rmspow",
    "ptpow",
    "nightpow",
    "timepow",
    "useorbMJD",
    "ecc_penalty",
    "max_astrom_rms",
    "maxrejnum",
)


class LinkMergeConfigTestCase(lsst.utils.tests.TestCase):

    def test_validate(self):
        LinkMergeConfig().validate()

    def test_shared_with_linkpurify(self):
        merge = LinkMergeConfig()
        purify = LinkPurifyConfig()
        for name in SHARED_KNOBS:
            self.assertEqual(
                getattr(merge, name),
                getattr(purify, name),
                msg=f"{name} drifted between LinkPurify and LinkMerge",
            )

    def test_max_oop_intentionally_tighter(self):
        # max_oop is intentionally tighter here (1000 vs purify's 10000).
        # Pinning so the divergence stays visible.
        self.assertEqual(LinkMergeConfig().max_oop, 1000.0)


class LinkMergeConnectionsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.connections = LinkMergeConnections(config=LinkMergeConfig())

    def test_dataset_type_names(self):
        all_conn = self.connections.allConnections
        self.assertEqual(all_conn["sspPurifiedLinkages"].name, "ssp_purified_linkages")
        self.assertEqual(all_conn["sspMergedLinkages"].name, "ssp_merged_linkages")
        self.assertEqual(
            all_conn["sspMergedLinkageSources"].name,
            "ssp_merged_linkage_sources",
        )

    def test_purified_inputs_multiple(self):
        # The point of this stage is to combine N balanced-index outputs
        # back into one table. The input needs multiple=True or the
        # quantum-graph builder only gathers one of them.
        self.assertTrue(self.connections.allConnections["sspPurifiedLinkages"].multiple)

    def test_quantum_dimensions(self):
        # Merge collapses the balanced-index axis.
        self.assertEqual(
            set(self.connections.dimensions),
            {"instrument", "day_obs", "ssp_hypothesis_table"},
        )


class LinkMergeTaskTestCase(lsst.utils.tests.TestCase):

    def test_default_name(self):
        self.assertEqual(LinkMergeTask._DefaultName, "linkMerge")

    def test_constructs(self):
        LinkMergeTask(config=LinkMergeConfig())

    def test_empty_input_raises_no_work_found(self):
        task = LinkMergeTask(config=LinkMergeConfig())
        with self.assertRaises(pipeBase.NoWorkFound):
            task.run(
                sspVisitInputs=None,
                sspTrackletSources=None,
                sspPurifiedLinkages=[],
                sspPurifiedLinkageSources=[],
            )


class LinkMergeRunTestCase(lsst.utils.tests.TestCase):
    """Run LinkMerge.run() against linkages produced by Heliolinc.

    Two copies in the list exercise the per-table clusternum renumber
    loop that sits in front of the heliolinx call.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        visits, sources, linkages, linkage_sources = make_linkage_inputs(bundle_id=1)
        cls.visits = visits
        cls.sources = sources
        cls.purified_linkages = [linkages, linkages.copy()]
        cls.purified_sources = [linkage_sources, linkage_sources.copy()]

    def test_run_pins_counts_and_schema(self):
        # Two copies of the same one-cluster table go in. The renumber loop
        # first bumps the second copy's clusternum to 1, then LinkMerge
        # deduplicates the overlap back to one cluster. Six sources stay
        # because only one copy's source rows map to the surviving cluster.
        result = LinkMergeTask(config=LinkMergeConfig()).run(
            self.visits, self.sources, self.purified_linkages, self.purified_sources,
        )
        self.assertEqual(len(result.sspMergedLinkages), 1)
        self.assertEqual(len(result.sspMergedLinkageSources), 6)
        self.assertEqual(set(result.sspMergedLinkages.dtype.names), EXPECTED_LINKAGE_COLS)
        self.assertEqual(int(result.sspMergedLinkages["clusternum"][0]), 0)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
