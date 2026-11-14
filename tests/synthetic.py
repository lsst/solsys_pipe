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

"""Small in-memory table builders used by the solsys_pipe tests.

Only the columns the production code actually reads are populated. That
way if someone renames a column upstream, the tests notice.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
from astropy.table import Table


# Arbitrary date inside the LSST survey window; not tied to a real visit.
DEFAULT_MJD = 60607.5

# Heliolinx 0.0.4 linkage columns. LinkPurify and LinkMerge only filter
# rows, so they hand back the same column set as Heliolinc itself.
EXPECTED_LINKAGE_COLS = {
    "astromRMS", "clusternum",
    "heliohyp0", "heliohyp1", "heliohyp2",
    "metric", "obsnights",
    "orbitVX", "orbitVY", "orbitVZ",
    "orbitX", "orbitY", "orbitZ",
    "orbit_MJD", "orbit_a", "orbit_e", "orbit_eval_count",
    "pairnum",
    "posRMS", "posX", "posY", "posZ",
    "rating", "reference_MJD", "timespan", "totRMS", "uniquepoints",
    "velRMS", "velX", "velY", "velZ",
}


def make_rng(seed=1905):
    """Return a deterministic numpy ``Generator``.

    Tests share one default seed so failures reproduce verbatim. Pass an
    explicit seed when a test needs a different draw without disturbing
    the shared fixtures.

    Parameters
    ----------
    seed : `int`, optional
        Seed for ``numpy.random.default_rng``.

    Returns
    -------
    rng : `numpy.random.Generator`
    """
    return np.random.default_rng(seed)


def make_dia_source_table(
    n_sources=8,
    n_visits=2,
    *,
    mjd0=DEFAULT_MJD,
    rng=None,
    with_reliability=True,
    with_trail=False,
):
    """Build a minimal DIA source table with pre-rename column names.

    Sources are split evenly across ``n_visits``; any remainder lands on
    the last visit. Visit MJDs are spaced one hour apart, well clear of
    the tracklet-pairing window heliolinx uses.

    Parameters
    ----------
    n_sources : `int`, optional
        Total number of source rows. Must be at least ``n_visits``.
    n_visits : `int`, optional
        Number of distinct visits.
    mjd0 : `float`, optional
        MJD of the first visit.
    rng : `numpy.random.Generator`, optional
        Defaults to ``make_rng()``.
    with_reliability : `bool`, optional
        Add a ``reliability`` column drawn uniformly on [0, 1].
    with_trail : `bool`, optional
        Add ``trailFlux``, ``trailLength``, and ``trailAngle`` columns.

    Returns
    -------
    table : `astropy.table.Table`
    """
    if rng is None:
        rng = make_rng()
    if n_visits < 1 or n_sources < n_visits:
        raise ValueError("need at least one source per visit")

    visit_ids = np.repeat(np.arange(1, n_visits + 1), n_sources // n_visits)
    if len(visit_ids) < n_sources:
        # Remainder lands on the last visit.
        visit_ids = np.concatenate(
            [visit_ids, np.full(n_sources - len(visit_ids), n_visits)]
        )
    visit_mjd = mjd0 + (visit_ids - 1) * (1.0 / 24.0)

    cols = {
        "visit": visit_ids,
        "diaSourceId": np.arange(1, n_sources + 1, dtype=np.int64),
        "midpointMjdTai": visit_mjd,
        "ra": rng.uniform(-0.01, 0.01, size=n_sources),
        "dec": rng.uniform(-0.01, 0.01, size=n_sources),
        "band": np.full(n_sources, "r", dtype="U1"),
        "psfFlux": rng.uniform(50.0, 500.0, size=n_sources),
        "psfFluxErr": rng.uniform(0.1, 1.0, size=n_sources),
        "raErr": np.full(n_sources, 1e-5),
        "decErr": np.full(n_sources, 1e-5),
    }
    if with_reliability:
        cols["reliability"] = rng.uniform(0.0, 1.0, size=n_sources)
    if with_trail:
        cols["trailFlux"] = rng.uniform(40.0, 400.0, size=n_sources)
        cols["trailLength"] = rng.uniform(0.0, 0.5, size=n_sources)
        cols["trailAngle"] = rng.uniform(0.0, 180.0, size=n_sources)
    return Table(cols)


def make_visit_summary_table(
    visit_ids: Iterable[int] | None = None,
    *,
    mjd0=DEFAULT_MJD,
    exposure_time=30.0,
):
    """Per-visit summary in the column form heliolinx expects.

    Parameters
    ----------
    visit_ids : iterable of `int`, optional
        Visit IDs; defaults to ``(1, 2)``.
    mjd0 : `float`, optional
        MJD of visit 1; later visits are spaced one hour apart.
    exposure_time : `float`, optional
        Seconds per exposure, identical for every visit.

    Returns
    -------
    table : `astropy.table.Table`
    """
    if visit_ids is None:
        visit_ids = (1, 2)
    visit_ids = np.asarray(list(visit_ids), dtype=np.int64)
    mjd = mjd0 + (visit_ids - 1) * (1.0 / 24.0)
    return Table({
        "visit": visit_ids,
        "MJD": mjd,
        "boresightRa": np.zeros(len(visit_ids)),
        "boresightDec": np.zeros(len(visit_ids)),
        "obsCode": np.full(len(visit_ids), "X05", dtype="U3"),
        "exposureTime": np.full(len(visit_ids), exposure_time),
    })


def make_pipeline_detections(
    n_nights=3,
    dets_per_night=2,
    *,
    mjd0=DEFAULT_MJD,
    ra0=10.0,
    dec0=-3.0,
    sky_rate_deg_per_day=0.05,
):
    """Build a DIA-style catalog of one slow-moving object.

    The returned tables match what ``MakeTrackletsTask.runQuantum`` hands
    to ``run()`` in production: columns are post-rename and carry the
    ``obscode`` we add to every detection.

    Sky rate stays well below heliolinx's default ``maxvel`` of
    1.5 deg/day, so tracklet pairs survive the velocity cut. Each
    detection sits in its own visit, with ``startind``/``endind`` pointing
    at the lone detection.

    Parameters
    ----------
    n_nights : `int`, optional
        Number of observation nights.
    dets_per_night : `int`, optional
        Detections per night.
    mjd0 : `float`, optional
        MJD of the first detection.
    ra0, dec0 : `float`, optional
        Sky position of the first detection in degrees.
    sky_rate_deg_per_day : `float`, optional
        Linear motion in RA.

    Returns
    -------
    detections : `astropy.table.Table`
        Per-detection table.
    visits : `astropy.table.Table`
        Per-visit table with Earth-state columns.
    """
    rows = []
    visit_rows = []
    i = 0
    for night in range(n_nights):
        for slot in range(dets_per_night):
            mjd = mjd0 + night + slot * (30.0 / 1440.0)  # 30 minutes apart.
            ra = ra0 + sky_rate_deg_per_day * (mjd - mjd0)
            rows.append({
                "MJD": mjd,
                "RA": ra,
                "Dec": dec0,
                "idstring": f"src_{i:03d}",
                "mag": 21.0,
                "band": "r",
                "trail_len": 0.0,
                "trail_PA": 0.0,
                "sig_across": 1e-4,
                "sig_along": 1e-4,
                "det_qual": 1000.0,
                "obscode": "X05",
            })
            # Earth pinned to (1, 0, 0) AU. Crude, but plausible enough
            # for heliolinx across this short observation window.
            visit_rows.append({
                "MJD": mjd,
                "RA": ra,
                "Dec": dec0,
                "obscode": "X05",
                "X": 1.0, "Y": 0.0, "Z": 0.0,
                "VX": 0.0, "VY": 0.0172, "VZ": 0.0,
                "startind": i,
                "endind": i + 1,
                "exptime": 30.0,
            })
            i += 1
    return Table(rows=rows), Table(rows=visit_rows)


def make_pipeline_hypothesis(bundle_id=1, r_au=3.0):
    """One-row hypothesis table for a main-belt-ish object.

    Parameters
    ----------
    bundle_id : `int`, optional
        SSP hypothesis bundle ID.
    r_au : `float`, optional
        Heliocentric radius in AU.

    Returns
    -------
    table : `astropy.table.Table`
    """
    return Table({
        "bundle_id": np.array([bundle_id], dtype=np.int64),
        "#r(AU)": np.array([r_au]),
        "rdot(AU/day)": np.array([0.0]),
        "mean_accel": np.array([0.0]),
    })


def make_pipeline_earth_state(mjd0=DEFAULT_MJD, n=10):
    """Tiny Earth ephemeris with the upstream uppercase column names.

    Parameters
    ----------
    mjd0 : `float`, optional
        MJD of the first sample.
    n : `int`, optional
        Number of samples, one per day.

    Returns
    -------
    table : `astropy.table.Table`
    """
    return Table({
        "MJD": mjd0 + np.arange(n, dtype=float),
        "X": np.full(n, 1.0),
        "Y": np.zeros(n),
        "Z": np.zeros(n),
        "VX": np.zeros(n),
        "VY": np.full(n, 0.0172),
        "VZ": np.zeros(n),
    })


def as_astropy_table(arr):
    """Wrap a numpy structured array as a Table (pass Tables through).

    Parameters
    ----------
    arr : `astropy.table.Table` or `numpy.ndarray`

    Returns
    -------
    table : `astropy.table.Table`
    """
    return arr if isinstance(arr, Table) else Table(arr)


def make_linkage_inputs(bundle_id=1):
    """Run MakeTracklets + Heliolinc on the synthetic detections.

    Heliolinc's output has the same dtype as the load-balancer's, so the
    linkage tables drop straight into LinkPurify or LinkMerge as their
    "balanced" or "purified" inputs.

    LSST imports happen inside the function so the module stays
    importable for quick ``python -c "import synthetic"`` smoke checks
    without the full stack.

    Parameters
    ----------
    bundle_id : `int`, optional
        SSP hypothesis bundle ID.

    Returns
    -------
    visits : `astropy.table.Table`
    sources : `astropy.table.Table`
    linkages : `astropy.table.Table`
    linkage_sources : `astropy.table.Table`
    """
    from lsst.solsys.pipe.heliolinc import HeliolincConfig, HeliolincTask
    from lsst.solsys.pipe.makeTracklets import MakeTrackletsConfig, MakeTrackletsTask

    detections, visits = make_pipeline_detections()
    sources, tracklets, t2s = MakeTrackletsTask(config=MakeTrackletsConfig()).run(
        detections, visits,
    )
    # MJDref must land inside the data window or heliolinx returns nothing.
    hel_config = HeliolincConfig()
    hel_config.MJDref = float(visits["MJD"].mean())
    result = HeliolincTask(config=hel_config).run(
        visits,
        as_astropy_table(sources),
        as_astropy_table(tracklets),
        as_astropy_table(t2s),
        make_pipeline_hypothesis(bundle_id=bundle_id),
        make_pipeline_earth_state(),
        ssp_hypothesis_bundle=bundle_id,
    )
    return (
        visits,
        as_astropy_table(sources),
        as_astropy_table(result.sspLinkage),
        as_astropy_table(result.sspLinkageSources),
    )


def make_tracklets_triplet(
    n_tracklets=3,
    sources_per_tracklet=2,
    *,
    visit_count=2,
    rng=None,
):
    """Return a (sources, tracklets, tracklet_to_source) triplet.

    Fills only the columns ConsolidateTracklets actually reads: index,
    image, Img1, Img2, trk_ID, i1, i2. The index mapping is dense (every
    tracklet and every source is referenced), which is the invariant
    ConsolidateTracklets asserts.

    Parameters
    ----------
    n_tracklets : `int`, optional
    sources_per_tracklet : `int`, optional
    visit_count : `int`, optional
        Distinct ``image`` values to spread the sources across.
    rng : `numpy.random.Generator`, optional
        Defaults to ``make_rng()``.

    Returns
    -------
    sources : `astropy.table.Table`
    tracklets : `astropy.table.Table`
    t2s : `astropy.table.Table`
        Tracklet-to-source index pairs.
    """
    if rng is None:
        rng = make_rng()
    n_sources = n_tracklets * sources_per_tracklet

    sources = Table({
        "index": np.arange(n_sources, dtype=np.int64),
        "image": np.repeat(np.arange(visit_count, dtype=np.int64),
                           max(1, n_sources // visit_count))[:n_sources],
        "MJD": rng.uniform(DEFAULT_MJD, DEFAULT_MJD + 0.5, size=n_sources),
        "RA": rng.uniform(-0.01, 0.01, size=n_sources),
        "Dec": rng.uniform(-0.01, 0.01, size=n_sources),
    })
    tracklets = Table({
        "trk_ID": np.arange(n_tracklets, dtype=np.int64),
        "Img1": np.zeros(n_tracklets, dtype=np.int64),
        "Img2": np.ones(n_tracklets, dtype=np.int64),
    })
    i1 = np.repeat(np.arange(n_tracklets, dtype=np.int64), sources_per_tracklet)
    i2 = np.arange(n_sources, dtype=np.int64)
    t2s = Table({"i1": i1, "i2": i2})
    return sources, tracklets, t2s
