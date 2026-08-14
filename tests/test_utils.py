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

"""Tests for lsst.solsys.pipe.utils."""

import unittest
from unittest import mock

import astropy.units as u
import numpy as np
from astropy.table import Table

import lsst.utils.tests
from lsst.solsys.pipe import utils

from synthetic import make_rng


# A throwaway dtype lets us exercise table_to_numpy without dragging the
# real heliolinx dtypes into the unit tests.
SIMPLE_DTYPE = np.dtype([("MJD", "f8"), ("RA", "f8"), ("Dec", "f8")])


def _random_simple_table(rng, size=8):
    """Random table with the columns of `SIMPLE_DTYPE`."""
    return Table({
        "MJD": rng.uniform(60000, 61000, size=size),
        "RA": rng.uniform(0, 360, size=size),
        "Dec": rng.uniform(-90, 90, size=size),
    })


class TableToNumpyTestCase(lsst.utils.tests.TestCase):

    def test_columns_round_trip(self):
        table = _random_simple_table(make_rng(seed=99))
        out = utils.table_to_numpy(table, SIMPLE_DTYPE)
        self.assertEqual(out.dtype, SIMPLE_DTYPE)
        self.assertEqual(len(out), len(table))
        for name in ("MJD", "RA", "Dec"):
            self.assertFloatsAlmostEqual(out[name], np.asarray(table[name]))

    def test_empty_table(self):
        # Quiet day_obs happens. A length-0 round trip has to work or the
        # NoWorkFound branch downstream is never reached.
        out = utils.table_to_numpy(Table({"MJD": [], "RA": [], "Dec": []}), SIMPLE_DTYPE)
        self.assertEqual(len(out), 0)
        self.assertEqual(out.dtype, SIMPLE_DTYPE)

    def test_missing_columns_default_to_zero(self):
        # When the source table doesn't carry a column the dtype expects,
        # the loop skips it and the field stays at its zero default. Pin
        # that so a future tightening doesn't surprise callers.
        partial = Table({"MJD": [60607.5, 60607.6], "RA": [10.0, 11.0]})
        out = utils.table_to_numpy(partial, SIMPLE_DTYPE)
        self.assertFloatsEqual(out["Dec"], np.zeros(2))

    def test_extra_columns_raise(self):
        # Flip side: a column not in the dtype raises ValueError. So
        # callers like MakeTracklets pre-trim before handing the table
        # to heliolinx.
        table = Table({
            "MJD": [60607.5], "RA": [10.0], "Dec": [-3.0],
            "psfFlux": [1.0],  # not in SIMPLE_DTYPE
        })
        with self.assertRaises(ValueError):
            utils.table_to_numpy(table, SIMPLE_DTYPE)


class RenameTableColumnsTestCase(lsst.utils.tests.TestCase):
    """Forgiving column-rename behavior."""

    def test_renames_present_columns_only(self):
        table = Table({"MJD": [1.0], "RA": [2.0], "Dec": [3.0]})
        utils.rename_table_columns(
            table,
            {
                "MJD": "mjd",          # present  -> renamed
                "missing": "ignored",  # absent   -> skipped
                "RA": "RA",            # identity -> skipped
            },
        )
        self.assertEqual(table.colnames, ["mjd", "RA", "Dec"])

    def test_empty_mapping_noop(self):
        table = Table({"MJD": [1.0], "RA": [2.0]})
        utils.rename_table_columns(table, {})
        self.assertEqual(table.colnames, ["MJD", "RA"])

    def test_colliding_targets_raise_and_leave_table_untouched(self):
        # Two columns both present and both mapped to "id" would end up
        # as duplicate columns. That must raise, name the offenders, and
        # leave the table exactly as it was.
        table = Table({"diaSourceId": [1], "sourceId": [2]})
        with self.assertRaisesRegex(ValueError, "diaSourceId.*sourceId.*'id'"):
            utils.rename_table_columns(
                table, {"diaSourceId": "id", "sourceId": "id"},
            )
        self.assertEqual(table.colnames, ["diaSourceId", "sourceId"])

    def test_collision_harmless_when_only_one_source_present(self):
        # The map may point several old names at the same new name, as
        # the DIA rename map does. That's fine as long as the table only
        # carries one of them.
        table = Table({"sourceId": [2], "RA": [3.0]})
        utils.rename_table_columns(
            table, {"diaSourceId": "id", "sourceId": "id"},
        )
        self.assertEqual(table.colnames, ["id", "RA"])


class TableColumnsToObjectArrayTestCase(lsst.utils.tests.TestCase):
    """Mixed-dtype columns survive the stack as Python objects."""

    def test_mixed_string_and_float_preserved(self):
        table = Table({
            "obsCode": ["I11", "I11"],
            "x": [1.5, 2.5],
            "y": [3.0, 4.0],
        })
        out = utils.table_columns_to_object_array(table, ["obsCode", "x", "y"])
        self.assertEqual(out.shape, (2, 3))
        self.assertEqual(out.dtype, object)
        self.assertEqual(out[0, 0], "I11")
        self.assertEqual(out[1, 1], 2.5)
        self.assertEqual(out[1, 2], 4.0)

    def test_quantity_columns_unwrap_to_floats(self):
        # The consolidateVisitTables=True path produces columns with
        # astropy units attached (boresightRa in degrees, etc.). Casting
        # through dtype=object drops the unit and leaves the bare numeric
        # value, which is what image_add_observerpos wants.
        table = Table({
            "MJD": [60607.5] * u.day,
            "boresightRa": [10.0] * u.deg,
            "obsCode": ["X05"],
        })
        out = utils.table_columns_to_object_array(table, ["MJD", "boresightRa", "obsCode"])
        self.assertEqual(out.shape, (1, 3))
        self.assertEqual(out[0, 0], 60607.5)
        self.assertEqual(out[0, 1], 10.0)
        self.assertEqual(out[0, 2], "X05")


class GroupedRangeMidpointsTestCase(lsst.utils.tests.TestCase):
    """Group midpoint math and output ordering."""

    def test_midpoints_and_group_order(self):
        # Rows are intentionally unsorted to confirm we group by key and
        # hand results back in ascending key order.
        table = Table({
            "visit": [20, 10, 20, 10, 30],
            "RA": [12.0, 5.0, 16.0, 9.0, -1.0],
            "Dec": [4.0, -3.0, 8.0, 1.0, 7.0],
        })
        ra_mid, dec_mid = utils.grouped_range_midpoints(
            table, "visit", ("RA", "Dec"),
        )
        # Output ordered by ascending group: 10, 20, 30.
        self.assertFloatsEqual(ra_mid, np.array([7.0, 14.0, -1.0]))
        self.assertFloatsEqual(dec_mid, np.array([-1.0, 6.0, 7.0]))

    def test_single_row_group(self):
        # A group with one row collapses to the value itself.
        table = Table({
            "visit": [1, 1, 2],
            "RA": [10.0, 14.0, 42.0],
            "Dec": [-1.0, 3.0, 5.0],
        })
        ra_mid, dec_mid = utils.grouped_range_midpoints(table, "visit", ("RA", "Dec"))
        self.assertFloatsEqual(ra_mid, np.array([12.0, 42.0]))
        self.assertFloatsEqual(dec_mid, np.array([1.0, 5.0]))


class TableToHeliolinxTestCase(lsst.utils.tests.TestCase):
    """Verify the factory lookup path and the column copy."""

    def test_uses_named_factory_and_copies_columns(self):
        table = Table({"MJD": [1.0, 2.0], "RA": [10.0, 11.0]})
        dtype = np.dtype([("MJD", "f8"), ("RA", "f8")])
        backing = np.zeros(2, dtype=dtype)

        # `create_demo` is a stand-in name; tell mock to add the attribute
        # for the duration of the patch rather than asserting it exists.
        with mock.patch.object(
            utils.hl, "create_demo", create=True, return_value=backing,
        ) as create_mock:
            out = utils.table_to_heliolinx(table, "demo")

        create_mock.assert_called_once_with(2)
        self.assertIs(out, backing)
        self.assertFloatsEqual(out["MJD"], np.asarray(table["MJD"]))
        self.assertFloatsEqual(out["RA"], np.asarray(table["RA"]))


class WrapperDelegationTestCase(lsst.utils.tests.TestCase):
    """make_hldet and make_hlimage are just thin wrappers."""

    def test_make_hldet_delegates_to_table_to_numpy(self):
        table = Table({"x": [1]})
        sentinel = object()
        with mock.patch.object(utils, "table_to_numpy", return_value=sentinel) as conv_mock:
            out = utils.make_hldet(table)
        self.assertIs(out, sentinel)
        conv_mock.assert_called_once_with(table, dtype=utils.solardg.hldet)

    def test_make_hlimage_delegates_to_table_to_numpy(self):
        table = Table({"x": [1]})
        sentinel = object()
        with mock.patch.object(utils, "table_to_numpy", return_value=sentinel) as conv_mock:
            out = utils.make_hlimage(table)
        self.assertIs(out, sentinel)
        conv_mock.assert_called_once_with(table, dtype=utils.solardg.hlimage)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
