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


import lsst.pipe.base as pipeBase
from astropy import units as u
from astropy.table import Table
from lsst.daf.base import DateTime
from lsst.pipe.tasks.postprocess import TableVStack


class ConsolidateSspTablesConnections(
    pipeBase.PipelineTaskConnections,
    defaultTemplates={
        "diaSourceInputName": "goodSeeingDiff_diaSrcTable",
        "diaSourceOutputName": "diaSourceRolling",
        "visitSummaryInputName": "finalVisitSummary",
        "visitInfoOutputName": "visitInfoRolling",
    },
    dimensions=("instrument", "day_obs"),
):
    inputDiaTables = pipeBase.connectionTypes.Input(
        doc="Input Source Tables to be concatenated",
        name="{diaSourceInputName}",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        multiple=True,
        deferLoad=True,
    )
    inputVisitSummaries = pipeBase.connectionTypes.Input(
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
        storageClass="ArrowAstropy",
        dimensions=("instrument", "day_obs"),
    )
    outputVisitInfo = pipeBase.connectionTypes.Output(
        doc="Concatenated Visit Summary from all day_obs in the input. The day_obs"
        "dimension would be the day_obs of the latest day_obs in the input.",
        name="{visitInfoOutputName}",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "day_obs"),
    )

    def adjust_all_quanta(self, adjuster):
        """This will drop all quanta but the quantum for the latest day_obs
        and add the input data from those quanta to the latest day_obs.
        """
        # Get all the data_ids to be iterated over.
        to_do = set(adjuster.iter_data_ids())

        # Dynamically get data_id for the latest day_obs.
        data_id_latest = max(to_do, key=lambda d: d["day_obs"])

        for data_id in to_do:
            # data_id has keys: instrument, day_obs.
            if data_id == data_id_latest:
                # Skip the latest day_obs. This is the one we want to keep.
                continue
            inputs = adjuster.get_inputs(data_id)
            # Since `multiple=True`, we need to loop over all input Refs.
            for input_data_id in inputs["inputDiaTables"]:
                # input_data_id has keys: instrument, visit, detector.
                # Add the input data taken from current data_id to the latest
                # day_obs.
                adjuster.add_input(data_id_latest, "inputDiaTables", input_data_id)
                # data_id_visit_summary has keys: instrument, visit.
                # We remove the detector key from data_id in the line below.
                data_id_visit_summary = input_data_id.subset(dimensions=self.inputVisitSummaries.dimensions)
                adjuster.add_input(data_id_latest, "inputVisitSummaries", data_id_visit_summary)
            adjuster.remove_quantum(data_id)


class ConsolidateSspTablesConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=ConsolidateSspTablesConnections
):
    pass


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

        assert len(inputDiaTableRefs) == len(inputVisitSummaryRefs), (
            f"Number of inputDiaTableRefs ({len(inputDiaTableRefs)}) does not match number of "
            f"inputVisitSummaryRefs ({len(inputVisitSummaryRefs)}). This is unexpected."
        )
        self.log.info(
            f"Concatenating {len(inputDiaTableRefs)} per-dayobs Source Tables "
            f"and {len(inputVisitSummaryRefs)} per-dayobs Visit Summaries"
        )
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

            # We ideally want the observatory code e.g. X05 for LSST, but
            # this is not available in visitInfo.

        # Make an Astropy table of visitInfo entries.
        consolidatedVisitInfo = Table(rows=ccdEntries)

        # Assign units to relevant columns.
        consolidatedVisitInfo["exposureTime"].unit = u.second
        consolidatedVisitInfo["MJD"].unit = u.day
        consolidatedVisitInfo["boresightRa"].unit = u.deg
        consolidatedVisitInfo["boresightDec"].unit = u.deg

        butlerQC.put(
            pipeBase.Struct(outputDiaTable=consolidatedDiaTable, outputVisitInfo=consolidatedVisitInfo),
            outputRefs,
        )
