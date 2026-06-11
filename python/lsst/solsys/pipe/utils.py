"""Helpers shared across solsys_pipe pipeline tasks."""
import heliolinx.heliolinx as hl
import heliolinx.solarsyst_dyn_geo as solardg
import numpy as np


def rename_table_columns(table, column_map):
    """Rename the columns that exist; quietly skip the rest.

    Astropy's ``rename_columns`` raises if a name in the map isn't on the
    table. This helper instead drops the missing names and renames the rest.
    We stay schema-forgiving because for now one rename map covers both the
    transient and DIA catalog schemas, and any given input only has some of
    those columns. That's safe because heliolinc can run on whatever survives.

    Parameters
    ----------
    table : `astropy.table.Table`
        Table to rename in place.
    column_map : `dict` [`str`, `str`]
        Maps old column name to new column name. Missing or identity
        entries are dropped silently.

    Raises
    ------
    ValueError
        If more than one column present on the table maps to the same new
        name, which would silently create duplicate columns.
    """
    renames = [(old, new) for old, new in column_map.items()
               if old != new and old in table.colnames]
    new_names = [new for _, new in renames]
    collisions = sorted({new for new in new_names if new_names.count(new) > 1})
    if collisions:
        offenders = {new: [old for old, mapped in renames if mapped == new]
                     for new in collisions}
        detail = "; ".join(
            f"{', '.join(old_cols)} would be renamed to '{new}'"
            for new, old_cols in offenders.items()
        )
        raise ValueError(
            f"Renaming would produce duplicate columns: {detail}. "
            f"The input should contain only one of them."
        )
    if renames:
        table.rename_columns(*zip(*renames))


def table_columns_to_object_array(table, column_names):
    """Stack columns side by side without losing per-column dtype.

    The heliolinx routine ``image_add_observerpos`` wants rows that mix a
    string observatory code with float coordinates. A plain
    ``np.column_stack`` picks one common dtype for the whole array and
    coerces floats into strings, which the C extension then rejects.
    Casting each column to ``object`` first sidesteps the coercion.

    Parameters
    ----------
    table : `astropy.table.Table`
        Source table.
    column_names : sequence of `str`
        Columns to stack, in output order.

    Returns
    -------
    array : `numpy.ndarray`
        Shape ``(len(table), len(column_names))`` with ``dtype=object``.
    """
    return np.column_stack(
        [np.asarray(table[name], dtype=object) for name in column_names]
    )


def grouped_range_midpoints(table, group_column, value_columns):
    """Return ``(min + max) / 2`` per group for each value column.

    Astropy's ``Table.group_by().groups.aggregate(...)`` would do the job
    but loops in Python once per group, which gets expensive once
    ``group_column`` has many thousands of distinct values (which it does
    on a real LSST night). We sort by the group key once and let
    ``np.minimum.reduceat`` / ``np.maximum.reduceat`` walk the sorted
    array in a single C pass.

    Parameters
    ----------
    table : `astropy.table.Table`
        Source table.
    group_column : `str`
        Column whose distinct values define the groups.
    value_columns : sequence of `str`
        Columns to reduce.

    Returns
    -------
    midpoints : `tuple` of `numpy.ndarray`
        One array per entry in ``value_columns``, each ordered by
        ascending group key.
    """
    keys = np.asarray(table[group_column])
    order = np.argsort(keys, kind="stable")
    # Where each new group starts in the sorted array.
    group_starts = np.r_[0, np.flatnonzero(np.diff(keys[order])) + 1]

    midpoints = []
    for name in value_columns:
        sorted_values = np.asarray(table[name], dtype=float)[order]
        lo = np.minimum.reduceat(sorted_values, group_starts)
        hi = np.maximum.reduceat(sorted_values, group_starts)
        midpoints.append(0.5 * (lo + hi))
    return tuple(midpoints)


def table_to_heliolinx(table, dtypename):
    """Copy a table into the heliolinx structured array named ``dtypename``.

    ``heliolinx.create_<dtypename>(n)`` hands back an empty structured
    array of the right dtype; we just fill each field from the matching
    table column.

    Parameters
    ----------
    table : `astropy.table.Table`
        Source table. Column names must match the heliolinx dtype fields.
    dtypename : `str`
        Suffix on the heliolinx factory function, e.g. ``"hlimage"``
        or ``"hldet"``.

    Returns
    -------
    array : `numpy.ndarray`
        Structured array of length ``len(table)``.
    """
    sa = getattr(hl, f"create_{dtypename}")(len(table))
    for name in table.colnames:
        sa[name] = np.asarray(table[name])
    return sa


def table_to_numpy(table, dtype):
    """Copy a table into a numpy structured array of the given dtype.

    Parameters
    ----------
    table : `astropy.table.Table`
        Source table.
    dtype : `numpy.dtype`
        Target structured dtype. Column names in ``table`` must be a
        subset of the dtype's fields; any field not present in the table
        stays at its zero default.

    Returns
    -------
    array : `numpy.ndarray`
        Structured array of length ``len(table)``.
    """
    sa = np.zeros(len(table), dtype=dtype)
    for name in table.colnames:
        sa[name] = np.asarray(table[name])
    return sa


def make_hldet(table):
    """Convert a table into the heliolinx ``hldet`` structured array.

    Parameters
    ----------
    table : `astropy.table.Table`
        Source table.

    Returns
    -------
    array : `numpy.ndarray`
        Structured array with dtype ``solardg.hldet``.
    """
    return table_to_numpy(table, dtype=solardg.hldet)


def make_hlimage(table):
    """Convert a table into the heliolinx ``hlimage`` structured array.

    Parameters
    ----------
    table : `astropy.table.Table`
        Source table.

    Returns
    -------
    array : `numpy.ndarray`
        Structured array with dtype ``solardg.hlimage``.
    """
    return table_to_numpy(table, dtype=solardg.hlimage)
