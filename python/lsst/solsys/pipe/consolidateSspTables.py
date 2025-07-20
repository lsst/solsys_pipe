# This file is part of solsys_pipe.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""A per-dayobs task which takes `goodSeeingDiff_assocDiaSrc` or
`goodSeeingDiff_diaSrcTable` from the last 14 (config-settable) days and
consolidates them into a `sspVisitInputs` as defined in `makeTracklets.py`
"""

__all__ = [
    "ConsolidateSspTablesTask",
    "ConsolidateSspTablesConfig",
]


import logging

import lsst.pex.config
import lsst.pipe.base as pipeBase
import numpy as np
from astropy import units as u
from astropy.table import Table
from heliolinx import solarsyst_dyn_geo as solardg
from lsst.daf.base import DateTime
from lsst.pipe.tasks.postprocess import TableVStack
from lsst.resources import ResourcePath
from lsst.solsys.pipe.utils import df2numpy

_LOG = logging.getLogger(__name__)


class ConsolidateSspTablesConnections(
    pipeBase.PipelineTaskConnections,
    defaultTemplates={
        "diaSourceInputName": "dia_source_detector",
        "earthStateInputName": "sspEarthState",
        "diaSourceOutputName": "dia_source_dayobs_14",
        "visitSummaryInputName": "visit_summary",
        "visitInfoOutputName": "visit_summary_dayobs_14",
    },
    dimensions=("instrument", "day_obs"),
):
    sspEarthState = pipeBase.connectionTypes.PrerequisiteInput(
        doc="Heliocentric Cartesian position and velocity for Earth",
        name="{earthStateInputName}",
        storageClass="ArrowAstropy",
        dimensions=(),
        deferLoad=True,
    )
    inputDiaTables = pipeBase.connectionTypes.PrerequisiteInput(
        doc="Input Source Tables to be concatenated",
        name="{diaSourceInputName}",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        multiple=True,
        deferLoad=True,
    )
    inputVisitSummaries = pipeBase.connectionTypes.PrerequisiteInput(
        doc="Per-visit consolidated exposure metadata",
        name="{visitSummaryInputName}",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
        multiple=True,
        deferLoad=True,
    )
    outputDiaTable = pipeBase.connectionTypes.Output(
        doc="Concatenated Source Table from all day_obs in the input. The day_obs"
        "dimension would be the day_obs of the latest day_obs in the input.",
        name="{diaSourceOutputName}",
        storageClass="DataFrame",
        dimensions=("instrument", "day_obs"),
    )
    outputVisitInfo = pipeBase.connectionTypes.Output(
        doc="Concatenated Visit Summary from all day_obs in the input. The day_obs"
        "dimension would be the day_obs of the latest day_obs in the input.",
        name="{visitInfoOutputName}",
        storageClass="DataFrame",
        dimensions=("instrument", "day_obs"),
    )

    def adjust_all_quanta(self, adjuster):
        """This will drop all quanta but the quantum for the latest day_obs
        and add the input data from those quanta to the latest day_obs.
        """
        # Get all the data_ids to be iterated over.
        to_do = set(adjuster.iter_data_ids())

        # If the iterable is empty, we have nothing to do.
        if not to_do:
            _LOG.warning("No data IDs to adjust quanta for.")
            return

        # Dynamically get data_id for the latest day_obs.
        data_id_latest = max(to_do, key=lambda d: d["day_obs"])

        # Loop over all data_id's in the to_do set. We will keep the latest
        # day_obs and add all the input data from the other day_obs to it.
        for data_id in to_do:
            # data_id has keys: instrument, day_obs.
            if data_id == data_id_latest:
                # Skip the latest day_obs. This is the one we want to keep.
                continue
            inputs = adjuster.get_prerequisite_inputs(data_id)
            # In the three for loops below, we loop over all input refs with
            # the same data_id as the current data_id in the loop and add
            # their input data to the quantum for the latest day_obs.
            for input_uuid in inputs["inputDiaTables"]:
                adjuster.add_prerequisite_input(data_id_latest, "inputDiaTables", input_uuid)
            for input_uuid in inputs["inputVisitSummaries"]:
                adjuster.add_prerequisite_input(data_id_latest, "inputVisitSummaries", input_uuid)
            adjuster.remove_quantum(data_id)

        # Log that the last day_obs is being kept as a reference.
        _LOG.info(
            f"Combined inputs from {len(to_do)} quanta into one quantum "
            f"under reference day_obs {data_id_latest['day_obs']}."
        )


class ConsolidateSspTablesConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=ConsolidateSspTablesConnections
):
    observatoryCode = lsst.pex.config.Field(
        dtype=str, default="X05", doc="Observatory code MPC, defaults to X05"
    )


class ConsolidateSspTablesTask(pipeBase.PipelineTask):
    """Concatenate `sourceTable` list into a per-dayobs `sourceTable_dayobs`"""

    ConfigClass = ConsolidateSspTablesConfig
    _DefaultName = "consolidateSspTables"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Concatenate the input source tables into a single output table."""
        inputs = butlerQC.get(inputRefs)
        inputDiaTableRefs = inputs["inputDiaTables"]
        inputVisitSummaryRefs = inputs["inputVisitSummaries"]

        # Let's make the order deterministic by sorting the input Refs.
        inputDiaTableRefs.sort(key=lambda x: (x.dataId["visit"], x.dataId["detector"]))

        # Sort by visit only since no detector key in visitSummary dataId.
        inputVisitSummaryRefs.sort(key=lambda x: x.dataId["visit"])

        # Visits in the visit summaries and unique visits in the DIA tables
        # should be the same.
        visitsInDia = set(ref.dataId["visit"] for ref in inputDiaTableRefs)
        visitsInSummary = set(ref.dataId["visit"] for ref in inputVisitSummaryRefs)
        # Some sanity checks.
        assert len(visitsInSummary) == len(inputVisitSummaryRefs), "Duplicate visits in visit summaries."
        assert visitsInDia == visitsInSummary, (
            f"Mismatch in visits: {len(visitsInSummary)} visits in visit summaries, "
            f"{len(visitsInDia)} unique visits in DIA tables."
        )

        self.log.info(
            f"Concatenating {len(inputVisitSummaryRefs)} visit summaries and {len(inputDiaTableRefs)} "
            f"DIA source tables over {len(visitsInDia)} visits.",
        )

        # Concatenate the input DIA tables into a single table without having
        # them all in memory at once.
        consolidatedDiaTable = TableVStack.vstack_handles(inputDiaTableRefs)

        # A list of dictionaries for each visit.
        ccdEntries = []

        for visitSummaryRef in inputVisitSummaryRefs:
            # visitInfo is the same for all elements in visitSummary.
            visitInfo = visitSummaryRef.get()[0].getVisitInfo()

            # Populate an entry with the visitInfo data.
            entry = dict(
                visit=visitInfo.id,  # == 'visit' in consolidatedDiaTable.
                exposureTime=visitInfo.exposureTime,
                MJD=visitInfo.date.get(system=DateTime.MJD),
                boresightRa=visitInfo.boresightRaDec[0].asDegrees(),
                boresightDec=visitInfo.boresightRaDec[1].asDegrees(),
                observatory=str(visitInfo.observatory),
                instrumentLabel=visitInfo.instrumentLabel,
                observationType=visitInfo.observationType,
                scienceProgram=visitInfo.scienceProgram,
                observationReason=visitInfo.observationReason,
                object=visitInfo.object,
                hasSimulatedContent=visitInfo.hasSimulatedContent,
            )
            ccdEntries.append(entry)

            # We ideally want the observatory code e.g. X05 for LSST, from
            # visitInfo, but since it's not available there, we retrieve it
            # from config instead (see below).

        # Make an Astropy table of visitInfo entries.
        consolidatedVisitInfo = Table(rows=ccdEntries)
        consolidatedVisitInfo["obsCode"] = self.config.observatoryCode
        image = (
            consolidatedVisitInfo[["MJD", "boresightRa", "boresightDec", "obsCode", "exposureTime"]]
            .to_pandas()
            .values
        )
        obsCodesTextLines = ResourcePath("resource://heliolinx/obsCodes.txt").read().decode().split("\n")
        obsarr = solardg.parse_ObsCodes(obsCodesTextLines)

        earthpos = df2numpy(
            inputs["sspEarthState"]
            .get()
            .to_pandas()
            .rename(columns={"X": "x", "Y": "y", "Z": "z", "VX": "vx", "VY": "vy", "VZ": "vz"}),
            "EarthState",
        )
        result = np.array(solardg.image_add_observerpos(image, obsarr, earthpos))

        for col in ["X", "Y", "Z", "VX", "VY", "VZ"]:
            consolidatedVisitInfo["observer" + col] = result[col]

        # Assign units to relevant columns.
        consolidatedVisitInfo["exposureTime"].unit = u.second
        consolidatedVisitInfo["MJD"].unit = u.day
        consolidatedVisitInfo["boresightRa"].unit = u.deg
        consolidatedVisitInfo["boresightDec"].unit = u.deg

        butlerQC.put(
            pipeBase.Struct(outputDiaTable=consolidatedDiaTable, outputVisitInfo=consolidatedVisitInfo),
            outputRefs,
        )
