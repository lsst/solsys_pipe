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

"""A per-dayobs task which consolidates `goodSeeingDiff_assocDiaSrc`,
`goodSeeingDiff_diaSrcTable` and consolidates them per-dayobs.`
"""

__all__ = [
    "ConsolidateSspDayObsTask",
    "ConsolidateSspDayObsConfig",
]


import logging

import lsst.pex.config
import lsst.pipe.base as pipeBase
import numpy as np
from . import utils
from astropy import units as u
from astropy.table import Table
from heliolinx import solarsyst_dyn_geo as solardg
from lsst.daf.base import DateTime
from lsst.pipe.tasks.postprocess import TableVStack
from lsst.resources import ResourcePath
from lsst.solsys.pipe.utils import df2numpy

_LOG = logging.getLogger(__name__)


class ConsolidateSspDayObsConnections(
    pipeBase.PipelineTaskConnections,
    defaultTemplates={
        "diaSourceInputName": "dia_source_visit",
        "earthStateInputName": "sspEarthState",
        "diaSourceOutputName": "dia_source_dayobs",
        "visitSummaryInputName": "visit_summary",
        "visitInfoOutputName": "visit_summary_dayobs",
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
        dimensions=("instrument", "visit"),
        multiple=True,
        minimum=0,
        deferLoad=True,
    )
    inputVisitSummaries = pipeBase.connectionTypes.PrerequisiteInput(
        doc="Per-visit consolidated exposure metadata",
        name="{visitSummaryInputName}",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
        multiple=True,
        minimum=0,
        deferLoad=True,
    )
    outputDiaTable = pipeBase.connectionTypes.Output(
        doc="Concatenated Source Table from one day_obs.",
        name="{diaSourceOutputName}",
        storageClass="DataFrame",
        dimensions=("instrument", "day_obs"),
    )
    outputVisitInfo = pipeBase.connectionTypes.Output(
        doc="Concatenated Visit Summary from one day_obs.",
        name="{visitInfoOutputName}",
        storageClass="DataFrame",
        dimensions=("instrument", "day_obs"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.consolidateVisitTables:
            self.prerequisiteInputs.remove("inputVisitSummaries")


diaSourceColumnRenameDict = {'diaSourceId': 'idstring', 'midpointMjdTai': 'MJD',
                             'ra': 'RA', 'dec': 'Dec', 'trailLength': 'trail_len', 'trailAngle': 'trail_PA',
                             'sourceId': 'idstring', 'expTime': 'exptime'}

class ConsolidateSspDayObsConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=ConsolidateSspDayObsConnections
):
    observatoryCode = lsst.pex.config.Field(
        dtype=str, default="X05", doc="Observatory code MPC, defaults to X05"
    )
    consolidateVisitTables = lsst.pex.config.Field(
        dtype=bool, default=True, doc="Whether to expect visit_summary-like data"
    )


class ConsolidateSspDayObsTask(pipeBase.PipelineTask):
    """Concatenate `sourceTable` list into a per-dayobs `sourceTable_dayobs`"""

    ConfigClass = ConsolidateSspDayObsConfig
    _DefaultName = "consolidateSspDayObs"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Concatenate the input source tables into a single output table."""
        inputs = butlerQC.get(inputRefs)
        inputDiaTableRefs = inputs["inputDiaTables"]
        if len(inputDiaTableRefs) == 0:
            raise pipeBase.NoWorkFound("No source tables found, skipping day_obs.")

        # Let's make the order deterministic by sorting the input Refs.
        inputDiaTableRefs.sort(key=lambda x: (x.dataId["visit"]))

        visitsInDia = set(ref.dataId["visit"] for ref in inputDiaTableRefs)

        self.log.info(
            f"Concatenating {len(inputDiaTableRefs)} DIA source tables over {len(visitsInDia)} visits.",
        )

        # Concatenate the input DIA tables into a single table without having
        # them all in memory at once.
        consolidatedDiaTable = TableVStack.vstack_handles(inputDiaTableRefs).to_pandas()
        consolidatedDiaTable = consolidatedDiaTable.rename(columns=diaSourceColumnRenameDict)
        consolidatedDiaTable = consolidatedDiaTable[['visit', 'MJD', 'RA', 'Dec', 'idstring']]
        consolidatedDiaTable['idstring'] = consolidatedDiaTable['idstring'].astype(str)
        consolidatedDiaTable['obscode'] = self.config.observatoryCode

        obsCodesTextLines = ResourcePath("resource://heliolinx/obsCodes.txt").read().decode().split("\n")
        obsarr = solardg.parse_ObsCodes(obsCodesTextLines)
        earthpos = df2numpy(
            inputs["sspEarthState"]
            .get()
            .to_pandas()
            .rename(columns={"X": "x", "Y": "y", "Z": "z", "VX": "vx", "VY": "vy", "VZ": "vz"}),
            "EarthState",
        )
        if self.config.consolidateVisitTables:
            inputVisitSummaryRefs = inputs["inputVisitSummaries"]
            inputVisitSummaryRefs.sort(key=lambda x: x.dataId["visit"])
            visitsInSummary = set(ref.dataId["visit"] for ref in inputVisitSummaryRefs)
            # Some sanity checks.
            assert len(visitsInSummary) == len(inputVisitSummaryRefs), "Duplicate visits in visit summaries."
            assert visitsInDia == visitsInSummary, (
                f"Mismatch in visits: {len(visitsInSummary)} visits in visit summaries, "
                f"{len(visitsInDia)} unique visits in DIA tables."
            )

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
            newimage = np.array(solardg.image_add_observerpos(image, obsarr, earthpos))

        else:
            def center(numbers):
                return (np.min(numbers) + np.max(numbers))/2
            groupby = consolidatedDiaTable[['visit', 'MJD', 'RA', 'Dec']].groupby('visit')
            mjd, ra, dec = groupby.aggregate(center).values.T
            mjd = mjd.astype(str).astype(float)
            ra = ra.astype(str).astype(float)
            dec = dec.astype(str).astype(float)
            obscode = np.repeat(self.config.observatoryCode, len(dec))
            expTime = np.repeat(30.0, len(dec))  # TODO: Make exact.
            image = np.array([mjd, ra, dec, obscode, expTime]).T
            newimage = solardg.image_add_observerpos(image, obsarr, earthpos)
        consolidatedVisitInfo = Table(newimage, names=['MJD', 'RA', 'Dec', 'obscode', 'X', 'Y', 'Z', 
                                                       'VX', 'VY', 'VZ', 'startind', 'endind', 'exptime'])
        # Assign units to relevant columns.
        consolidatedVisitInfo["exptime"].unit = u.second
        consolidatedVisitInfo["MJD"].unit = u.day
        consolidatedVisitInfo["RA"].unit = u.deg
        consolidatedVisitInfo["Dec"].unit = u.deg

        consolidatedDiaTable = consolidatedDiaTable[['MJD', 'RA', 'Dec', 'idstring']]

        butlerQC.put(
            pipeBase.Struct(outputDiaTable=consolidatedDiaTable, outputVisitInfo=consolidatedVisitInfo),
            outputRefs,
        )
