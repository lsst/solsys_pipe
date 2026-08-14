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

"""Tests for ConsolidateTracklets.

Two things matter here. ``run()`` stitches per-dayobs tables together
and has to keep all the cross-table indices consistent - that's plain
astropy work driven by synthetic tables. ``adjust_all_quanta`` decides
which day_obs quanta share a linking window, so we stand up a real
butler with the SSP dimensions and let AllDimensionsQuantumGraphBuilder
run the same way pipetask would (the pattern is borrowed from
``pipe_base/tests/test_adjust_all_quanta.py``).
"""

import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
from astropy.table import Table

import lsst.daf.butler.tests as butlerTests
import lsst.pipe.base as pipeBase
from lsst.pipe.base import PipelineGraph
from lsst.pipe.base.all_dimensions_quantum_graph_builder import (
    AllDimensionsQuantumGraphBuilder,
)
import lsst.utils.tests

from lsst.solsys.pipe.consolidateTracklets import (
    ConsolidateTrackletsConfig,
    ConsolidateTrackletsConnections,
    ConsolidateTrackletsTask,
)

from butler_setup import make_ssp_repo, populate_ssp_dimensions
from synthetic import make_tracklets_triplet, make_rng


# These tuples mirror what ConsolidateTrackletsConnections declares:
# visit_summary_* uses (instrument, day_obs); the tracklet datasets add
# ssp_hypothesis_table. The QG builder rejects any mismatch, so when
# the connections grow or drop a dimension, update here too.
_INPUT_DATASET_TYPES = (
    ("visit_summary_dayobs",
     {"instrument", "day_obs"}, "ArrowAstropy"),
    ("ssp_tracklet_source_dayobs",
     {"instrument", "day_obs", "ssp_hypothesis_table"}, "ArrowAstropy"),
    ("ssp_tracklet_dayobs",
     {"instrument", "day_obs", "ssp_hypothesis_table"}, "ArrowAstropy"),
    ("ssp_tracklet_to_source_dayobs",
     {"instrument", "day_obs", "ssp_hypothesis_table"}, "ArrowAstropy"),
)
_OUTPUT_DATASET_TYPES = (
    ("visit_summary_dayobs_14",
     {"instrument", "day_obs"}, "ArrowAstropy"),
    ("ssp_tracklet_source_dayobs_14",
     {"instrument", "day_obs", "ssp_hypothesis_table"}, "ArrowAstropy"),
    ("ssp_tracklet_dayobs_14",
     {"instrument", "day_obs", "ssp_hypothesis_table"}, "ArrowAstropy"),
    ("ssp_tracklet_to_source_dayobs_14",
     {"instrument", "day_obs", "ssp_hypothesis_table"}, "ArrowAstropy"),
)


class ConsolidateConfigTestCase(lsst.utils.tests.TestCase):

    def test_defaults(self):
        # 14-day window is the canonical SSP value, baked into the
        # pipeline YAMLs and the dataset-type suffixes.
        self.assertEqual(ConsolidateTrackletsConfig().linkingTimespan, 14)

    def test_validate(self):
        ConsolidateTrackletsConfig().validate()


class ConsolidateTaskTestCase(lsst.utils.tests.TestCase):

    def test_default_name(self):
        self.assertEqual(ConsolidateTrackletsTask._DefaultName, "consolidateTracklets")

    def test_constructs(self):
        ConsolidateTrackletsTask(config=ConsolidateTrackletsConfig())


class ConsolidateConnectionsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.connections = ConsolidateTrackletsConnections(
            config=ConsolidateTrackletsConfig()
        )

    def test_dataset_type_names(self):
        # These names are what connect daily tracklets to the 14-day
        # tracklet tables used by Heliolinc.
        all_conn = self.connections.allConnections
        self.assertEqual(all_conn["inputVisitSummaries"].name,
                         "visit_summary_dayobs")
        self.assertEqual(all_conn["inputTrackletSources"].name,
                         "ssp_tracklet_source_dayobs")
        self.assertEqual(all_conn["inputTracklets"].name,
                         "ssp_tracklet_dayobs")
        self.assertEqual(all_conn["inputTrackletToSource"].name,
                         "ssp_tracklet_to_source_dayobs")
        self.assertEqual(all_conn["outputVisitSummaries"].name,
                         "visit_summary_dayobs_14")
        self.assertEqual(all_conn["outputTrackletSources"].name,
                         "ssp_tracklet_source_dayobs_14")
        self.assertEqual(all_conn["outputTracklets"].name,
                         "ssp_tracklet_dayobs_14")
        self.assertEqual(all_conn["outputTrackletToSource"].name,
                         "ssp_tracklet_to_source_dayobs_14")

    def test_quantum_dimensions(self):
        self.assertEqual(
            set(self.connections.dimensions),
            {"instrument", "day_obs", "ssp_hypothesis_table"},
        )


class ConsolidateRunQuantumTestCase(lsst.utils.tests.TestCase):
    """Check runQuantum sorts day_obs refs before it reads inputs."""

    def _ref(self, day_obs):
        ref = mock.Mock()
        ref.dataId = {"day_obs": day_obs}
        return ref

    def test_sorts_all_inputs_by_day_obs_before_get(self):
        task = ConsolidateTrackletsTask(config=ConsolidateTrackletsConfig())
        inputRefs = mock.Mock()
        for name in ("inputVisitSummaries", "inputTrackletSources",
                     "inputTracklets", "inputTrackletToSource"):
            setattr(inputRefs, name, [self._ref(20240103), self._ref(20240101)])

        butlerQC = mock.Mock()
        butlerQC.get.return_value = {
            "inputVisitSummaries": [],
            "inputTrackletSources": [],
            "inputTracklets": [],
            "inputTrackletToSource": [],
        }
        outputs = object()

        with mock.patch.object(task, "run", return_value=outputs) as run_mock:
            task.runQuantum(butlerQC, inputRefs, outputRefs="outputs")

        # The important bit here is the handoff order: Butler get sees
        # sorted refs, then those fetched inputs go straight into run().
        butlerQC.get.assert_called_once_with(inputRefs)
        run_mock.assert_called_once_with(**butlerQC.get.return_value)
        butlerQC.put.assert_called_once_with(outputs, "outputs")
        for name in ("inputVisitSummaries", "inputTrackletSources",
                     "inputTracklets", "inputTrackletToSource"):
            self.assertEqual(
                [ref.dataId["day_obs"] for ref in getattr(inputRefs, name)],
                [20240101, 20240103],
            )


class DayobsToMjdTestCase(lsst.utils.tests.TestCase):
    """Tiny helper, used inside the window picker."""

    def setUp(self):
        self.connections = ConsolidateTrackletsConnections(
            config=ConsolidateTrackletsConfig()
        )

    def test_known_date(self):
        # 2024-01-01 is MJD 60310.
        self.assertEqual(self.connections.dayobs_to_mjd("20240101"), 60310)

    def test_accepts_int(self):
        # Some butler builds hand day_obs in as an int - the cast survives.
        self.assertEqual(self.connections.dayobs_to_mjd(20240101), 60310)

    def test_consecutive_days_increment_by_one(self):
        # Two days back to back must give MJDs that differ by one. Any
        # calendar drift in the conversion shows up here.
        a = self.connections.dayobs_to_mjd("20240101")
        b = self.connections.dayobs_to_mjd("20240102")
        self.assertEqual(b - a, 1)

    def test_year_boundary(self):
        # 2023-12-31 -> 2024-01-01 also one MJD apart; checks the
        # year rollover.
        a = self.connections.dayobs_to_mjd("20231231")
        b = self.connections.dayobs_to_mjd("20240101")
        self.assertEqual(b - a, 1)

    def test_leap_day(self):
        # 2024 is a leap year. Feb 28 -> Feb 29 -> Mar 1 each one MJD
        # apart. Catches anyone tempted to "simplify" with a 30-days-
        # per-month shortcut.
        feb28 = self.connections.dayobs_to_mjd("20240228")
        feb29 = self.connections.dayobs_to_mjd("20240229")
        mar01 = self.connections.dayobs_to_mjd("20240301")
        self.assertEqual(feb29 - feb28, 1)
        self.assertEqual(mar01 - feb29, 1)

    def test_wrong_length_raises(self):
        with self.assertRaises(AssertionError):
            self.connections.dayobs_to_mjd("2024-01")


class ConsolidateRunTestCase(lsst.utils.tests.TestCase):
    """Index renumber checks for the consolidator."""

    def setUp(self):
        self.task = ConsolidateTrackletsTask(config=ConsolidateTrackletsConfig())
        # Three days, different sizes, so cumulative offsets aren't trivial.
        self.day_inputs = [
            make_tracklets_triplet(
                n_tracklets=2 + i,
                sources_per_tracklet=2,
                rng=make_rng(seed=2000 + i),
            )
            for i in range(3)
        ]
        # Visit summaries pass through without changes. A trivial stub works.
        self.visit_summaries = [
            Table({
                "MJD": np.array([60000.0 + i, 60000.1 + i]),
                "visit": np.array([2 * i + 1, 2 * i + 2], dtype=np.int64),
            })
            for i in range(3)
        ]

    def _split(self):
        return (
            [t[0] for t in self.day_inputs],
            [t[1] for t in self.day_inputs],
            [t[2] for t in self.day_inputs],
        )

    def test_total_lengths(self):
        sources, tracklets, t2s = self._split()
        res = self.task.run(self.visit_summaries, sources, tracklets, t2s)
        self.assertEqual(len(res.outputTracklets), sum(len(t) for t in tracklets))
        self.assertEqual(len(res.outputTrackletSources), sum(len(s) for s in sources))
        self.assertEqual(len(res.outputTrackletToSource), sum(len(t) for t in t2s))
        self.assertEqual(len(res.outputVisitSummaries),
                         sum(len(v) for v in self.visit_summaries))

    def test_trk_ids_dense(self):
        # Stacked trk_IDs should run 0..N-1 with no gaps.
        sources, tracklets, t2s = self._split()
        res = self.task.run(self.visit_summaries, sources, tracklets, t2s)
        ids = np.asarray(res.outputTracklets["trk_ID"])
        self.assertFloatsEqual(ids, np.arange(len(ids)))

    def test_t2s_indices_in_bounds(self):
        # i1/i2 must land inside the concatenated tracklet/source tables.
        sources, tracklets, t2s = self._split()
        res = self.task.run(self.visit_summaries, sources, tracklets, t2s)
        self.assertLess(res.outputTrackletToSource["i1"].max(),
                        len(res.outputTracklets))
        self.assertLess(res.outputTrackletToSource["i2"].max(),
                        len(res.outputTrackletSources))
        self.assertGreaterEqual(res.outputTrackletToSource["i1"].min(), 0)
        self.assertGreaterEqual(res.outputTrackletToSource["i2"].min(), 0)

    def test_reproducible(self):
        # Same input, same output. Easy to miss, easy to check.
        sources, tracklets, t2s = self._split()
        a = self.task.run(
            [v.copy() for v in self.visit_summaries],
            [s.copy() for s in sources],
            [t.copy() for t in tracklets],
            [t.copy() for t in t2s],
        )
        b = self.task.run(
            [v.copy() for v in self.visit_summaries],
            [s.copy() for s in sources],
            [t.copy() for t in tracklets],
            [t.copy() for t in t2s],
        )
        self.assertFloatsEqual(a.outputTracklets["trk_ID"],
                               b.outputTracklets["trk_ID"])
        self.assertFloatsEqual(a.outputTrackletToSource["i1"],
                               b.outputTrackletToSource["i1"])
        self.assertFloatsEqual(a.outputTrackletToSource["i2"],
                               b.outputTrackletToSource["i2"])

    def test_single_day_passes_through(self):
        # One day in, one day out. Total counts must match exactly and
        # trk_IDs must already be 0..N-1, since there's no cross-day
        # renumbering to do.
        s, t, t2 = make_tracklets_triplet(
            n_tracklets=3, sources_per_tracklet=2, rng=make_rng(seed=1),
        )
        visit_one = [self.visit_summaries[0]]
        res = self.task.run(visit_one, [s], [t], [t2])
        self.assertEqual(len(res.outputTracklets), len(t))
        self.assertEqual(len(res.outputTrackletSources), len(s))
        self.assertEqual(len(res.outputTrackletToSource), len(t2))
        ids = np.asarray(res.outputTracklets["trk_ID"])
        self.assertFloatsEqual(ids, np.arange(len(t)))

    def test_second_day_indices_offset_by_first_day_size(self):
        # Two-day input. The second day's t2s indices should land at
        # offsets equal to the first day's sizes, since the consolidator
        # concatenates and renumbers in order.
        s0, t0, t2_0 = make_tracklets_triplet(
            n_tracklets=2, sources_per_tracklet=2, rng=make_rng(seed=10),
        )
        s1, t1, t2_1 = make_tracklets_triplet(
            n_tracklets=3, sources_per_tracklet=2, rng=make_rng(seed=11),
        )
        res = self.task.run(self.visit_summaries[:2], [s0, s1], [t0, t1],
                            [t2_0, t2_1])
        # First day's t2s i1 stays in [0, len(t0)). Second day's i1 should
        # all be >= len(t0).
        i1 = np.asarray(res.outputTrackletToSource["i1"])
        i2 = np.asarray(res.outputTrackletToSource["i2"])
        self.assertGreaterEqual(int(i1[len(t2_0):].min()), len(t0))
        self.assertLess(int(i1[:len(t2_0)].max()), len(t0))
        self.assertGreaterEqual(int(i2[len(t2_0):].min()), len(s0))
        self.assertLess(int(i2[:len(t2_0)].max()), len(s0))

    def test_empty_raises_no_work_found(self):
        with self.assertRaises(pipeBase.NoWorkFound):
            self.task.run([], [], [], [])

    def test_mismatched_lengths_caught(self):
        sources, tracklets, t2s = self._split()
        with self.assertRaises(AssertionError):
            self.task.run(self.visit_summaries[:-1], sources, tracklets, t2s)


class AdjustAllQuantaTestCase(lsst.utils.tests.TestCase):
    """Drive adjust_all_quanta the way pipetask does.

    Repo gets built once per class. Doing it per test costs a few
    seconds we don't need. Each test gets its own collection so butler
    puts don't leak across tests.
    """

    INSTRUMENT = "LSSTCam"
    HYPOTHESIS_TABLE = "default"
    # Four days. With linkingTimespan=14 the first three roll up into
    # one window, the fourth sits ~29 days later so it's on its own.
    DAY_OBS = (20240101, 20240102, 20240103, 20240201)

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.root = tempfile.mkdtemp()
        cls.repo = make_ssp_repo(cls.root)
        populate_ssp_dimensions(
            cls.repo,
            instrument=cls.INSTRUMENT,
            day_obs_values=cls.DAY_OBS,
            hypothesis_table=cls.HYPOTHESIS_TABLE,
        )
        for name, dims, storage_class in _INPUT_DATASET_TYPES:
            butlerTests.addDatasetType(cls.repo, name, dims, storage_class)
        for name, dims, storage_class in _OUTPUT_DATASET_TYPES:
            butlerTests.addDatasetType(cls.repo, name, dims, storage_class)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)
        super().tearDownClass()

    def setUp(self):
        super().setUp()
        self.butler = butlerTests.makeTestCollection(self.repo, uniqueId=self.id())
        # Seed empty tables at every (day_obs, hypothesis_table) so the
        # quantum-graph builder has something to plan against. Content
        # doesn't matter, adjust_all_quanta only looks at data IDs.
        empty_visit = Table({"MJD": np.array([], dtype=float)})
        empty_tracklet_like = Table({"trk_ID": np.array([], dtype=np.int64)})
        for day in self.DAY_OBS:
            base_id = {"instrument": self.INSTRUMENT, "day_obs": day}
            ssp_id = dict(base_id, ssp_hypothesis_table=self.HYPOTHESIS_TABLE)
            self.butler.put(empty_visit, "visit_summary_dayobs", base_id)
            self.butler.put(empty_tracklet_like, "ssp_tracklet_source_dayobs", ssp_id)
            self.butler.put(empty_tracklet_like, "ssp_tracklet_dayobs", ssp_id)
            self.butler.put(
                Table({"i1": np.array([], dtype=np.int64),
                       "i2": np.array([], dtype=np.int64)}),
                "ssp_tracklet_to_source_dayobs", ssp_id,
            )

    def _build_quanta(self, linking_timespan):
        """Build a QG and return {day_obs: Quantum} for the survivors."""
        config = ConsolidateTrackletsConfig()
        config.linkingTimespan = linking_timespan

        graph = PipelineGraph(universe=self.butler.dimensions)
        graph.add_task("consolidateTracklets", ConsolidateTrackletsTask, config=config)

        qgb = AllDimensionsQuantumGraphBuilder(
            graph,
            self.butler,
            input_collections=[self.butler.run],
            output_run=self.butler.run + "_out",
        )
        qg = qgb.finish(attach_datastore_records=False).assemble()
        quanta = qg.build_execution_quanta(task_label="consolidateTracklets")
        return {q.dataId["day_obs"]: q for q in quanta.values()}

    def test_consecutive_days_collapse_into_latest(self):
        # 14-day window: the three Jan days roll up into Jan 3 and Feb 1
        # sits on its own.
        quanta = self._build_quanta(linking_timespan=14)
        self.assertEqual(set(quanta), {20240103, 20240201})

        latest = quanta[20240103]
        # Four multi-input connections across three Jan days = 12 refs
        # flowing into the survivor quantum.
        total_input_refs = sum(len(refs) for refs in latest.inputs.values())
        self.assertGreaterEqual(total_input_refs, 12)

    def test_short_window_keeps_disjoint_days(self):
        # 1-day window: only adjacent days link. Jan 1-3 still collapse
        # to Jan 3 because they're consecutive; Feb 1 is on its own.
        quanta = self._build_quanta(linking_timespan=1)
        self.assertIn(20240103, quanta)
        self.assertIn(20240201, quanta)

    def test_quantum_dimensions(self):
        quanta = self._build_quanta(linking_timespan=14)
        for q in quanta.values():
            self.assertEqual(
                set(q.dataId.dimensions.names),
                {"instrument", "day_obs", "ssp_hypothesis_table"},
            )


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
