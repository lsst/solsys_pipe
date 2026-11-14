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

"""Helpers for tests that need a real butler with the SSP dimensions.

solsys_pipe defines its own dimensions (``ssp_hypothesis_table``,
``ssp_hypothesis_bundle``, ``ssp_balanced_index``) in
``prototypes/dimensions.yaml``. We hand that file to ``makeTestRepo`` so
the test butler knows about them. Same idea as the existing daf_butler
tests, just with a custom ``DimensionConfig``.
"""

from __future__ import annotations

from contextlib import contextmanager
import os
import warnings

from astropy.utils.metadata import MergeConflictWarning
from lsst.daf.butler import DimensionConfig
import lsst.daf.butler.tests as butlerTests


# Walk up from this file to the package root, then into prototypes/.
_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(_TESTS_DIR)
DIMENSIONS_YAML = os.path.join(PACKAGE_DIR, "prototypes", "dimensions.yaml")


def ssp_dimension_config():
    """Standard butler universe with the SSP-specific dimensions added.

    Returns
    -------
    config : `lsst.daf.butler.DimensionConfig`
    """
    return DimensionConfig(DIMENSIONS_YAML)


def make_ssp_repo(root):
    """Create a fresh in-memory test repo that knows the SSP dimensions.

    Parameters
    ----------
    root : `str`
        Directory for the test repo's on-disk artifacts.

    Returns
    -------
    repo : `lsst.daf.butler.Butler`
    """
    return butlerTests.makeTestRepo(root, dimensionConfig=ssp_dimension_config())


@contextmanager
def ignore_butler_metadata_merge_warnings():
    """Silence the harmless ``LSST.BUTLER.*`` merge-conflict warnings.

    The test Butler stamps every input table with its own
    ``LSST.BUTLER.*`` provenance (dataset IDs, data IDs). When a task
    stacks input tables, astropy warns because the combined table can
    only keep one value per metadata key. The pipeline tests assert on
    rows and columns, not on which input's provenance survives, so the
    warnings are pure noise.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Cannot merge meta key 'LSST\.BUTLER\..*",
            category=MergeConflictWarning,
        )
        yield


def populate_ssp_dimensions(
    repo,
    *,
    instrument="LSSTCam",
    day_obs_values=(20240101, 20240102, 20240103),
    hypothesis_table="default",
    hypothesis_bundle_ids=(),
    balanced_index_ids=(),
):
    """Seed the dimension records SSP tasks need to plan quanta against.

    Parameters
    ----------
    repo : `lsst.daf.butler.Butler`
        Test repo produced by `make_ssp_repo`.
    instrument : `str`, optional
        Instrument label registered with the butler.
    day_obs_values : iterable of `int`, optional
        ``day_obs`` values to register under ``instrument``.
    hypothesis_table : `str`, optional
        SSP hypothesis table label.
    hypothesis_bundle_ids : iterable of `int`, optional
        Bundle IDs to register against ``hypothesis_table``.
    balanced_index_ids : iterable of `int`, optional
        Balanced-index IDs to register against ``hypothesis_table``.
    """
    butlerTests.addDataIdValue(repo, "instrument", instrument)
    butlerTests.addDataIdValue(repo, "ssp_hypothesis_table", hypothesis_table)
    for day in day_obs_values:
        butlerTests.addDataIdValue(repo, "day_obs", day, instrument=instrument)
    for bundle in hypothesis_bundle_ids:
        butlerTests.addDataIdValue(
            repo, "ssp_hypothesis_bundle", bundle,
            ssp_hypothesis_table=hypothesis_table,
        )
    for idx in balanced_index_ids:
        butlerTests.addDataIdValue(
            repo, "ssp_balanced_index", idx,
            ssp_hypothesis_table=hypothesis_table,
        )
