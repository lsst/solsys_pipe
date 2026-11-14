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

"""Construction, config, and connection checks for MakeTracklets.

The real tracklet finding lives inside heliolinx. What we cover from
the Python side is the glue: config defaults, the DIA-to-heliolinx
rename map, and the connection schema downstream YAMLs depend on.
"""

import unittest
from unittest import mock

import astropy.units as u
import numpy as np
from astropy.table import Table
import lsst.pipe.base as pipeBase
from lsst.pipe.base import InMemoryDatasetHandle
import lsst.utils.tests
from lsst.solsys.pipe import utils
from lsst.solsys.pipe.makeTracklets import (
    MakeTrackletsConfig,
    MakeTrackletsConnections,
    MakeTrackletsTask,
    diaSourceColumnRenameDict,
)

from synthetic import make_dia_source_table, make_pipeline_detections, make_pipeline_earth_state


class MakeTrackletsConfigTestCase(lsst.utils.tests.TestCase):

    def test_defaults(self):
        # Pin survey/observatory facts and structural switches only.
        # Science tuning knobs (maxvel, matchrad, maxgcr, ...) belong to
        # the science team - pinning them here means every retune is a
        # spurious test failure.
        c = MakeTrackletsConfig()
        self.assertEqual(c.observatoryCode, "X05")
        self.assertEqual(c.exptime, 30)
        self.assertEqual(c.imagerad, 2.0)
        self.assertEqual(c.mintrkpts, 2)
        self.assertTrue(c.consolidateVisitTables)

    def test_validate(self):
        MakeTrackletsConfig().validate()


class MakeTrackletsRenameMapTestCase(lsst.utils.tests.TestCase):
    """The contract between DIA columns and the heliolinx names."""

    def test_known_columns_route_correctly(self):
        self.assertEqual(diaSourceColumnRenameDict["midpointMjdTai"], "MJD")
        self.assertEqual(diaSourceColumnRenameDict["ra"], "RA")
        self.assertEqual(diaSourceColumnRenameDict["dec"], "Dec")
        self.assertEqual(diaSourceColumnRenameDict["diaSourceId"], "idstring")
        self.assertEqual(diaSourceColumnRenameDict["sourceId"], "idstring")

    def test_target_set_complete(self):
        # Check nothing was dropped from the map.
        expected = {
            "idstring", "MJD", "RA", "Dec",
            "trail_len", "trail_PA",
            "exptime", "sig_along", "sig_across",
        }
        self.assertEqual(set(diaSourceColumnRenameDict.values()), expected)


class MakeTrackletsAstropyTablePrepTestCase(lsst.utils.tests.TestCase):
    """Small checks for the Astropy-table prep done before heliolinx."""

    def test_dia_column_rename_stays_astropy(self):
        table = make_dia_source_table(n_sources=4, n_visits=2)
        utils.rename_table_columns(table, diaSourceColumnRenameDict)

        self.assertIsInstance(table, Table)
        for name in ("idstring", "MJD", "RA", "Dec", "sig_along", "sig_across"):
            self.assertIn(name, table.colnames)
        for name in ("diaSourceId", "midpointMjdTai", "ra", "dec"):
            self.assertNotIn(name, table.colnames)

    def test_duplicate_target_column_names_are_rejected(self):
        table = Table({
            "diaSourceId": [1, 2],
            "sourceId": [10, 20],
        })
        with self.assertRaises(ValueError):
            utils.rename_table_columns(table, diaSourceColumnRenameDict)

    def test_visit_centers_match_old_grouping_contract(self):
        # Visits come in unsorted, but the grouped output is sorted by visit.
        table = Table({
            "visit": [2, 1, 2, 1],
            "MJD": [10.0, 1.0, 14.0, 5.0],
            "RA": [100.0, 10.0, 120.0, 20.0],
            "Dec": [0.0, -5.0, 10.0, 5.0],
        })
        mjd, ra, dec = utils.grouped_range_midpoints(
            table,
            "visit",
            ("MJD", "RA", "Dec"),
        )

        self.assertFloatsEqual(mjd, np.array([3.0, 12.0]))
        self.assertFloatsEqual(ra, np.array([15.0, 110.0]))
        self.assertFloatsEqual(dec, np.array([0.0, 5.0]))

    def test_mixed_columns_become_rows_for_observer_position_code(self):
        table = Table({
            "MJD": [60607.5],
            "RA": [10.0],
            "obsCode": ["X05"],
        })
        rows = utils.table_columns_to_object_array(table, ["MJD", "RA", "obsCode"])

        self.assertEqual(rows.shape, (1, 3))
        self.assertEqual(rows[0, 2], "X05")


class MakeTrackletsTaskTestCase(lsst.utils.tests.TestCase):

    def test_default_name(self):
        self.assertEqual(MakeTrackletsTask._DefaultName, "makeTracklets")

    def test_constructs(self):
        MakeTrackletsTask(config=MakeTrackletsConfig())

    def test_no_visit_summary_drops_input(self):
        # When consolidateVisitTables is False the pipeline doesn't have
        # per-visit summaries to pass in, so the connection has to drop
        # off the inputs list. Otherwise the quantum-graph builder
        # complains that an input is missing.
        c = MakeTrackletsConfig()
        c.consolidateVisitTables = False
        connections = MakeTrackletsConnections(config=c)
        self.assertNotIn("inputVisitSummaries", connections.inputs)


class MakeTrackletsConnectionsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.connections = MakeTrackletsConnections(config=MakeTrackletsConfig())

    def test_dataset_type_names(self):
        # These names are the contract other pipeline steps use.
        all_conn = self.connections.allConnections
        self.assertEqual(all_conn["sspEarthState"].name, "sspEarthState")
        self.assertEqual(all_conn["inputDiaTables"].name, "dia_source_visit")
        self.assertEqual(all_conn["inputVisitSummaries"].name, "visit_summary")
        self.assertEqual(all_conn["outputDiaTable"].name, "dia_source_dayobs")
        self.assertEqual(all_conn["outputDiaTable"].storageClass, "ArrowAstropy")
        self.assertEqual(all_conn["outputVisitInfo"].name, "visit_summary_dayobs")
        self.assertEqual(all_conn["outputVisitInfo"].storageClass, "ArrowAstropy")
        self.assertEqual(
            all_conn["sspTrackletSources"].name,
            "ssp_tracklet_source_dayobs",
        )
        self.assertEqual(all_conn["sspTracklets"].name, "ssp_tracklet_dayobs")
        self.assertEqual(
            all_conn["sspTrackletToSource"].name,
            "ssp_tracklet_to_source_dayobs",
        )

    def test_inputs_and_outputs(self):
        all_conn = self.connections.allConnections
        for name in ("sspEarthState", "inputDiaTables", "inputVisitSummaries"):
            self.assertIn(name, all_conn)
        for name in ("outputDiaTable", "outputVisitInfo",
                     "sspTrackletSources", "sspTracklets", "sspTrackletToSource"):
            self.assertIn(name, all_conn)

    def test_quantum_dimensions(self):
        # SSP runs per instrument, day_obs and hypothesis table.
        self.assertEqual(
            set(self.connections.dimensions),
            {"instrument", "day_obs", "ssp_hypothesis_table"},
        )


class MakeTrackletsRunQuantumTestCase(lsst.utils.tests.TestCase):
    """Check runQuantum's Butler-facing table prep before heliolinx."""

    def _dia_table(self, visit, source_ids, mjd, ra, dec, psf_flux, reliability):
        n = len(source_ids)
        return Table({
            "visit": np.full(n, visit, dtype=np.int64),
            "diaSourceId": np.asarray(source_ids, dtype=np.int64),
            "midpointMjdTai": np.asarray(mjd, dtype=float),
            "ra": np.asarray(ra, dtype=float),
            "dec": np.asarray(dec, dtype=float),
            "band": np.full(n, "r", dtype="U1"),
            "psfFlux": np.asarray(psf_flux, dtype=float),
            "psfFluxErr": np.full(n, 1.0),
            "raErr": np.full(n, 1.0e-5),
            "decErr": np.full(n, 1.0e-5),
            "reliability": np.asarray(reliability, dtype=float),
        })

    @staticmethod
    def _fake_image_add_observerpos(image, obsarr, earthpos):
        # Preserve the first five columns from runQuantum's image array;
        # fill observer state columns with simple dummy values.
        rows = []
        for i, row in enumerate(image):
            rows.append([
                float(row[0]), float(row[1]), float(row[2]), str(row[3]),
                1.0, 2.0, 3.0, 0.1, 0.2, 0.3, i, i + 1, float(row[4]),
            ])
        return np.asarray(rows, dtype=object)

    def _run_quantum(self, task, butler_inputs):
        """Run runQuantum with the heliolinx boundary mocked out.

        Returns the tables handed to run() plus the mocks, so tests can
        check the prep work without invoking the C extension.
        """
        butlerQC = mock.Mock()
        butlerQC.get.return_value = butler_inputs

        prepared_tables = {}
        fake_outputs = (
            Table({"index": np.array([42], dtype=np.int64)}),
            Table({"trk_ID": np.array([43], dtype=np.int64)}),
            Table({
                "i1": np.array([0], dtype=np.int64),
                "i2": np.array([0], dtype=np.int64),
            }),
        )

        def record_run(inputDiaTables, inputVisitSummaries):
            # Stop before heliolinx. These tests care about the tables
            # runQuantum prepared and handed to run().
            prepared_tables["dia"] = inputDiaTables
            prepared_tables["visits"] = inputVisitSummaries
            return fake_outputs

        obs_codes_resource = mock.Mock()
        obs_codes_resource.read.return_value = b""

        with mock.patch("lsst.solsys.pipe.makeTracklets.ResourcePath",
                        return_value=obs_codes_resource), \
                mock.patch("lsst.solsys.pipe.makeTracklets.solardg.parse_ObsCodes",
                           return_value=object()), \
                mock.patch("lsst.solsys.pipe.makeTracklets.solardg.image_add_observerpos",
                           side_effect=self._fake_image_add_observerpos), \
                mock.patch("lsst.solsys.pipe.makeTracklets.utils.table_to_heliolinx",
                           return_value=object()), \
                mock.patch.object(task, "run", side_effect=record_run) as run_mock:
            task.runQuantum(butlerQC, inputRefs=None, outputRefs="outputs")

        return prepared_tables, butlerQC, run_mock, fake_outputs

    def test_run_quantum_prepares_astropy_tables_without_visit_summaries(self):
        config = MakeTrackletsConfig()
        config.consolidateVisitTables = False
        config.minReliability = 0.5
        good_reliability = (1.0 + config.minReliability) / 2.0
        bad_reliability = config.minReliability / 2.0
        bad_flux = -1.0
        task = MakeTrackletsTask(config=config)

        # Both visit-2 rows pass the cuts; their midpoint lands at
        # MJD=12, RA=110, Dec=5 in outputVisitInfo.
        dia_visit_2 = self._dia_table(
            2,
            [3, 4],
            [10.0, 14.0],
            [100.0, 120.0],
            [0.0, 10.0],
            [100.0, 200.0],
            [good_reliability, good_reliability],
        )
        # First two visit-1 rows pass (midpoint: MJD=3, RA=15, Dec=0).
        # Row 3 is cut by reliability; row 4 by negative flux.
        dia_visit_1 = self._dia_table(
            1,
            [1, 2, 5, 6],
            [1.0, 5.0, 999.0, 777.0],
            [10.0, 20.0, 999.0, 777.0],
            [-5.0, 5.0, 999.0, 777.0],
            [100.0, 200.0, 300.0, bad_flux],
            [good_reliability, good_reliability, bad_reliability, good_reliability],
        )
        prepared_tables, butlerQC, run_mock, fake_outputs = self._run_quantum(task, {
            # Deliberately unsorted, so this test checks runQuantum sorts
            # deferred DIA inputs before stacking them.
            "inputDiaTables": [
                InMemoryDatasetHandle(dia_visit_2, storageClass="ArrowAstropy", visit=2),
                InMemoryDatasetHandle(dia_visit_1, storageClass="ArrowAstropy", visit=1),
            ],
            "sspEarthState": InMemoryDatasetHandle(
                make_pipeline_earth_state(n=2),
                storageClass="ArrowAstropy",
            ),
        })
        fake_sources, fake_tracklets, fake_t2s = fake_outputs

        butlerQC.get.assert_called_once_with(None)
        run_mock.assert_called_once()
        dia_table = prepared_tables["dia"]
        visit_table = prepared_tables["visits"]

        self.assertIsInstance(dia_table, Table)
        self.assertIsInstance(visit_table, Table)
        self.assertEqual(len(dia_table), 4)
        self.assertNotIn("visit", dia_table.colnames)
        for name in ("idstring", "MJD", "RA", "Dec", "mag", "obscode", "det_qual"):
            self.assertIn(name, dia_table.colnames)
        for name in ("diaSourceId", "midpointMjdTai", "ra", "dec"):
            self.assertNotIn(name, dia_table.colnames)
        self.assertEqual(list(np.asarray(dia_table["idstring"], dtype=str)), ["1", "2", "3", "4"])

        self.assertFloatsEqual(np.asarray(visit_table["MJD"], dtype=float), np.array([3.0, 12.0]))
        self.assertFloatsEqual(np.asarray(visit_table["RA"], dtype=float), np.array([15.0, 110.0]))
        self.assertFloatsEqual(np.asarray(visit_table["Dec"], dtype=float), np.array([0.0, 5.0]))
        self.assertFloatsEqual(np.asarray(visit_table["exptime"], dtype=float), np.array([30.0, 30.0]))
        self.assertEqual(list(np.asarray(visit_table["obscode"], dtype=str)), ["X05", "X05"])
        butlerQC.put.assert_called_once()
        outputs, output_refs = butlerQC.put.call_args.args
        self.assertIs(outputs.outputDiaTable, dia_table)
        self.assertIs(outputs.outputVisitInfo, visit_table)
        self.assertIs(outputs.sspTrackletSources, fake_sources)
        self.assertIs(outputs.sspTracklets, fake_tracklets)
        self.assertIs(outputs.sspTrackletToSource, fake_t2s)
        self.assertEqual(output_refs, "outputs")

    def test_no_dia_tables_raises_no_work_found(self):
        # A day_obs with no DIA tables should skip cleanly, not crash.
        task = MakeTrackletsTask(config=MakeTrackletsConfig())
        butlerQC = mock.Mock()
        butlerQC.get.return_value = {"inputDiaTables": []}
        with self.assertRaises(pipeBase.NoWorkFound):
            task.runQuantum(butlerQC, inputRefs=None, outputRefs=None)

    def test_trail_flux_preferred_with_psf_fallback(self):
        # mag comes from trailFlux when the trail fit worked, and falls
        # back to psfFlux where the fit returned NaN. The positive-flux
        # cut applies to whichever flux was chosen.
        config = MakeTrackletsConfig()
        config.consolidateVisitTables = False
        task = MakeTrackletsTask(config=config)

        dia = self._dia_table(
            1,
            [1, 2, 3],
            [1.0, 1.0, 1.0],
            [10.0, 11.0, 12.0],
            [0.0, 0.0, 0.0],
            [100.0, 50.0, 300.0],
            [1.0, 1.0, 1.0],
        )
        # Row 1: failed trail fit -> psfFlux=100. Row 2: good trail
        # fit -> trailFlux=200. Row 3: negative trailFlux -> dropped even
        # though its psfFlux is fine.
        dia["trailFlux"] = np.array([np.nan, 200.0, -5.0])
        dia["trailLength"] = np.zeros(3)
        dia["trailAngle"] = np.zeros(3)

        prepared_tables, _, _, _ = self._run_quantum(task, {
            "inputDiaTables": [
                InMemoryDatasetHandle(dia, storageClass="ArrowAstropy", visit=1),
            ],
            "sspEarthState": InMemoryDatasetHandle(
                make_pipeline_earth_state(n=2),
                storageClass="ArrowAstropy",
            ),
        })

        dia_table = prepared_tables["dia"]
        self.assertEqual(len(dia_table), 2)
        expected_mag = np.array([
            (100.0 * u.nJy).to_value(u.ABmag),
            (200.0 * u.nJy).to_value(u.ABmag),
        ])
        self.assertFloatsAlmostEqual(
            np.asarray(dia_table["mag"], dtype=float), expected_mag, rtol=1e-12,
        )

    def _visit_summary_ref(self, visit, mjd, ra, dec, exptime=30.0):
        """Mock one deferred visit_summary handle with just enough
        visitInfo for runQuantum to read."""
        visitInfo = mock.Mock()
        visitInfo.id = visit
        visitInfo.exposureTime = exptime
        visitInfo.date.get.return_value = mjd
        visitInfo.boresightRaDec = (
            mock.Mock(asDegrees=mock.Mock(return_value=ra)),
            mock.Mock(asDegrees=mock.Mock(return_value=dec)),
        )
        visitInfo.observatory = "LSST"
        visitInfo.instrumentLabel = "LSSTCam"
        visitInfo.observationType = "science"
        visitInfo.scienceProgram = "survey"
        visitInfo.observationReason = "science"
        visitInfo.object = ""
        visitInfo.hasSimulatedContent = False
        ref = mock.Mock()
        ref.dataId = {"visit": visit}
        ref.get.return_value = [mock.Mock(getVisitInfo=mock.Mock(return_value=visitInfo))]
        return ref

    def test_run_quantum_builds_visit_table_from_visit_summaries(self):
        # With consolidateVisitTables=True (the default) the visit table
        # comes from each visit's visitInfo, not from DIA midpoints.
        task = MakeTrackletsTask(config=MakeTrackletsConfig())
        dia_visit_1 = self._dia_table(
            1, [1, 2], [1.0, 5.0], [10.0, 20.0], [-5.0, 5.0],
            [100.0, 200.0], [1.0, 1.0],
        )
        dia_visit_2 = self._dia_table(
            2, [3, 4], [10.0, 14.0], [100.0, 120.0], [0.0, 10.0],
            [100.0, 200.0], [1.0, 1.0],
        )
        prepared_tables, _, run_mock, _ = self._run_quantum(task, {
            "inputDiaTables": [
                InMemoryDatasetHandle(dia_visit_1, storageClass="ArrowAstropy", visit=1),
                InMemoryDatasetHandle(dia_visit_2, storageClass="ArrowAstropy", visit=2),
            ],
            # Unsorted on purpose: runQuantum sorts summaries by visit.
            "inputVisitSummaries": [
                self._visit_summary_ref(2, mjd=12.5, ra=110.0, dec=5.0),
                self._visit_summary_ref(1, mjd=3.5, ra=15.0, dec=0.0),
            ],
            "sspEarthState": InMemoryDatasetHandle(
                make_pipeline_earth_state(n=2),
                storageClass="ArrowAstropy",
            ),
        })

        run_mock.assert_called_once()
        visit_table = prepared_tables["visits"]
        self.assertFloatsEqual(np.asarray(visit_table["MJD"], dtype=float),
                               np.array([3.5, 12.5]))
        self.assertFloatsEqual(np.asarray(visit_table["RA"], dtype=float),
                               np.array([15.0, 110.0]))
        self.assertFloatsEqual(np.asarray(visit_table["Dec"], dtype=float),
                               np.array([0.0, 5.0]))
        self.assertFloatsEqual(np.asarray(visit_table["exptime"], dtype=float),
                               np.array([30.0, 30.0]))
        self.assertEqual(list(np.asarray(visit_table["obscode"], dtype=str)),
                         ["X05", "X05"])


class MakeTrackletsRunTestCase(lsst.utils.tests.TestCase):
    """Actually run MakeTracklets.run() and look at what comes back."""

    def setUp(self):
        # 3 nights, 2 detections each, ~30 min apart. Sky rate slow.
        self.detections, self.visits = make_pipeline_detections()
        self.task = MakeTrackletsTask(config=MakeTrackletsConfig())

    # Schema heliolinx 0.0.4 hands back. Pinned from a real run so a
    # field rename or addition trips the test instead of silently
    # showing up downstream.
    EXPECTED_SOURCE_COLS = {
        "Dec",
        "MJD",
        "RA",
        "band",
        "det_qual",
        "idstring",
        "image",
        "index",
        "known_obj",
        "mag",
        "obscode",
        "sig_across",
        "sig_along",
        "sigmag",
        "trail_PA",
        "trail_len",
    }
    EXPECTED_TRACKLET_COLS = {
        "Dec1",
        "Dec2",
        "Img1",
        "Img2",
        "RA1",
        "RA2",
        "npts",
        "trk_ID",
    }
    EXPECTED_T2S_COLS = {"i1", "i2"}

    def test_finds_tracklets(self):
        # 3 nights with 2 detections each, close in time and slow on the
        # sky. Heliolinx pairs each night into one tracklet, so we get
        # 6 sources, 3 tracklets, 6 t2s rows (two members per tracklet).
        sources, tracklets, t2s = self.task.run(self.detections, self.visits)
        self.assertEqual(len(sources), 6)
        self.assertEqual(len(tracklets), 3)
        self.assertEqual(len(t2s), 6)

    def test_output_schemas(self):
        # Column sets are the API between heliolinx and the rest of SSP.
        # If they ever drift, we want to know before downstream tasks
        # silently miss a field.
        sources, tracklets, t2s = self.task.run(self.detections, self.visits)
        self.assertEqual(set(sources.dtype.names), self.EXPECTED_SOURCE_COLS)
        self.assertEqual(set(tracklets.dtype.names), self.EXPECTED_TRACKLET_COLS)
        self.assertEqual(set(t2s.dtype.names), self.EXPECTED_T2S_COLS)

    def test_t2s_indices_cover_full_range(self):
        # Every tracklet should have its members in t2s, and every source
        # should be referenced. So i1 spans 0..n_tracklets-1 and i2 spans
        # 0..n_sources-1.
        sources, tracklets, t2s = self.task.run(self.detections, self.visits)
        self.assertEqual(int(t2s["i1"].min()), 0)
        self.assertEqual(int(t2s["i1"].max()), len(tracklets) - 1)
        self.assertEqual(int(t2s["i2"].min()), 0)
        self.assertEqual(int(t2s["i2"].max()), len(sources) - 1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
