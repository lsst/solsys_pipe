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

"""Config and connection checks for LinkPurify."""

import unittest

import lsst.pipe.base as pipeBase
import lsst.utils.tests
from lsst.solsys.pipe.linkPurify import (
    LinkPurifyConfig,
    LinkPurifyConnections,
    LinkPurifyTask,
)

from synthetic import EXPECTED_LINKAGE_COLS, make_linkage_inputs


class LinkPurifyConfigTestCase(lsst.utils.tests.TestCase):

    def test_validate(self):
        LinkPurifyConfig().validate()

    def test_planarity_backend_default_on(self):
        # Flipping doLinkPlanarity changes which heliolinx function gets
        # called, so it's structural. The shared tuning knobs themselves
        # are checked over in test_linkMerge.py.
        self.assertTrue(LinkPurifyConfig().doLinkPlanarity)


class LinkPurifyConnectionsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.connections = LinkPurifyConnections(config=LinkPurifyConfig())

    def test_dataset_type_names(self):
        all_conn = self.connections.allConnections
        self.assertEqual(all_conn["sspBalancedLinkages"].name, "ssp_balanced_linkages")
        self.assertEqual(all_conn["sspBalancedLinkageSources"].name,
                         "ssp_balanced_linkage_sources")
        self.assertEqual(all_conn["sspPurifiedLinkages"].name, "ssp_purified_linkages")
        self.assertEqual(all_conn["sspPurifiedLinkageSources"].name,
                         "ssp_purified_linkage_sources")

    def test_quantum_dimensions(self):
        # LinkPurify runs per balanced-index slot from the load balancer.
        self.assertEqual(
            set(self.connections.dimensions),
            {"instrument", "day_obs", "ssp_hypothesis_table", "ssp_balanced_index"},
        )


class LinkPurifyTaskTestCase(lsst.utils.tests.TestCase):

    def test_default_name(self):
        self.assertEqual(LinkPurifyTask._DefaultName, "linkPurify")

    def test_constructs(self):
        LinkPurifyTask(config=LinkPurifyConfig())

    def test_empty_input_raises_no_work_found(self):
        # The empty-input branch has to short-circuit before reaching
        # heliolinx, otherwise the C extension throws on a zero-row input.
        task = LinkPurifyTask(config=LinkPurifyConfig())
        with self.assertRaises(pipeBase.NoWorkFound):
            task.run(
                sspVisitInputs=None,
                sspTrackletSources=None,
                sspBalancedLinkages=[],
                sspBalancedLinkageSources=[],
            )


class LinkPurifyRunTestCase(lsst.utils.tests.TestCase):
    """Run LinkPurify.run() against linkages produced by Heliolinc.

    Heliolinc's output has the same column layout as the load balancer's, so
    these linkages can go straight into LinkPurify without running the balancer.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.visits, cls.sources, cls.balanced_linkages, cls.balanced_sources = (
            make_linkage_inputs(bundle_id=1)
        )

    def test_run_pins_counts_and_schema(self):
        # Single clean cluster passes the quality filters unchanged.
        result = LinkPurifyTask(config=LinkPurifyConfig()).run(
            self.visits, self.sources, self.balanced_linkages, self.balanced_sources,
        )
        self.assertEqual(len(result.sspPurifiedLinkages), 1)
        self.assertEqual(len(result.sspPurifiedLinkageSources), 6)
        self.assertEqual(set(result.sspPurifiedLinkages.dtype.names), EXPECTED_LINKAGE_COLS)
        self.assertEqual(int(result.sspPurifiedLinkages["clusternum"][0]), 0)

    def test_classic_linkpurify_backend_agrees(self):
        # doLinkPlanarity=False switches to the older linkPurify routine.
        # A single clean cluster should survive either backend the same way.
        config = LinkPurifyConfig()
        config.doLinkPlanarity = False
        result = LinkPurifyTask(config=config).run(
            self.visits, self.sources, self.balanced_linkages, self.balanced_sources,
        )
        self.assertEqual(len(result.sspPurifiedLinkages), 1)
        self.assertEqual(len(result.sspPurifiedLinkageSources), 6)
        self.assertEqual(set(result.sspPurifiedLinkages.dtype.names), EXPECTED_LINKAGE_COLS)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
