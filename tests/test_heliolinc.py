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

"""Config and connection checks for HeliolincTask.

The linker itself lives in heliolinx. What we cover here is the LSST
side: defaults, dataset names, dimensions.
"""

import unittest
from types import SimpleNamespace
from unittest import mock

import lsst.utils.tests
from lsst.solsys.pipe.heliolinc import (
    HeliolincConfig,
    HeliolincConnections,
    HeliolincTask,
)
from lsst.solsys.pipe.makeTracklets import MakeTrackletsConfig, MakeTrackletsTask

from synthetic import (
    EXPECTED_LINKAGE_COLS,
    as_astropy_table,
    make_pipeline_detections,
    make_pipeline_earth_state,
    make_pipeline_hypothesis,
)


class HeliolincConfigTestCase(lsst.utils.tests.TestCase):

    def test_validate(self):
        HeliolincConfig().validate()

    def test_geocentric_grid_well_ordered(self):
        # Sanity check: max has to be bigger than min or the grid search
        # has nothing to walk through.
        c = HeliolincConfig()
        self.assertGreater(c.maxgeodist, c.mingeodist)

    def test_bound_orbits_by_default(self):
        # max_v_inf == 0 keeps the search to bound orbits only. Flipping
        # this changes what objects can be found, so it's a real choice
        # worth pinning.
        c = HeliolincConfig()
        self.assertEqual(c.use_univar, 0)
        self.assertEqual(c.max_v_inf, 0.0)


class HeliolincConnectionsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.connections = HeliolincConnections(config=HeliolincConfig())

    def test_dataset_type_names(self):
        all_conn = self.connections.allConnections
        self.assertEqual(all_conn["sspTracklets"].name, "ssp_tracklet_dayobs_14")
        self.assertEqual(all_conn["sspTrackletSources"].name,
                         "ssp_tracklet_source_dayobs_14")
        self.assertEqual(all_conn["sspTrackletToSource"].name,
                         "ssp_tracklet_to_source_dayobs_14")
        self.assertEqual(all_conn["sspLinkage"].name, "ssp_linkages")
        self.assertEqual(all_conn["sspLinkageSources"].name, "ssp_linkage_sources")

    def test_quantum_dimensions(self):
        self.assertEqual(
            set(self.connections.dimensions),
            {"instrument", "day_obs", "ssp_hypothesis_table", "ssp_hypothesis_bundle"},
        )


class HeliolincTaskTestCase(lsst.utils.tests.TestCase):

    def test_default_name(self):
        self.assertEqual(HeliolincTask._DefaultName, "link")

    def test_constructs(self):
        HeliolincTask(config=HeliolincConfig())


class HeliolincRunQuantumTestCase(lsst.utils.tests.TestCase):
    """Check runQuantum passes the quantum's bundle id into run()."""

    def test_passes_hypothesis_bundle_id_to_run(self):
        task = HeliolincTask(config=HeliolincConfig())
        butlerQC = mock.Mock()
        butlerQC.quantum = SimpleNamespace(
            dataId=SimpleNamespace(
                ssp_hypothesis_bundle=SimpleNamespace(id=7),
            ),
        )
        butlerQC.get.return_value = {"sspVisitInputs": object()}
        outputs = object()

        with mock.patch.object(task, "run", return_value=outputs) as run_mock:
            task.runQuantum(butlerQC, inputRefs="inputs", outputRefs="outputs")

        # The only custom logic here is adding ssp_hypothesis_bundle=7.
        butlerQC.get.assert_called_once_with("inputs")
        run_mock.assert_called_once_with(
            sspVisitInputs=butlerQC.get.return_value["sspVisitInputs"],
            ssp_hypothesis_bundle=7,
        )
        butlerQC.put.assert_called_once_with(outputs, "outputs")


class HeliolincRunTestCase(lsst.utils.tests.TestCase):
    """Run Heliolinc.run() with the tracklet output from MakeTracklets."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        detections, visits = make_pipeline_detections()
        mt = MakeTrackletsTask(config=MakeTrackletsConfig())
        sources, tracklets, t2s = mt.run(detections, visits)
        # Wrap numpy outputs as Tables to match what the butler hands downstream.
        cls.visits = visits
        cls.sources = as_astropy_table(sources)
        cls.tracklets = as_astropy_table(tracklets)
        cls.t2s = as_astropy_table(t2s)
        cls.hypothesis = make_pipeline_hypothesis(bundle_id=1)
        cls.earth = make_pipeline_earth_state()

    def _make_task(self):
        # heliolinx's MJDref defaults to 0 (year 1858), well outside our data
        # window. With the wrong reference time it returns no linkages and no
        # error, so tests would silently pass on empty output. Use the midpoint.
        config = HeliolincConfig()
        config.MJDref = float(self.visits["MJD"].mean())
        return HeliolincTask(config=config)

    def test_run_pins_counts_and_schema(self):
        # One synthetic object observed over 3 nights produces one cluster of
        # 3 tracklets. Each of the 6 detections gets its own source-membership
        # row in sspLinkageSources.
        result = self._make_task().run(
            self.visits, self.sources, self.tracklets, self.t2s,
            self.hypothesis, self.earth, ssp_hypothesis_bundle=1,
        )
        self.assertEqual(len(result.sspLinkage), 1)
        self.assertEqual(len(result.sspLinkageSources), 6)
        self.assertEqual(set(result.sspLinkage.dtype.names), EXPECTED_LINKAGE_COLS)
        self.assertEqual(set(result.sspLinkageSources.dtype.names), {"i1", "i2"})
        self.assertEqual(int(result.sspLinkage["clusternum"][0]), 0)
        self.assertEqual(int(result.sspLinkageSources["i1"].min()), 0)
        self.assertEqual(int(result.sspLinkageSources["i1"].max()), 0)
        self.assertEqual(int(result.sspLinkageSources["i2"].min()), 0)
        self.assertEqual(int(result.sspLinkageSources["i2"].max()), 5)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
