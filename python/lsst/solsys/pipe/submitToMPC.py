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


from heliolinx import solarsyst_dyn_geo as solardg
import logging
import numpy as np

import lsst.pipe.base as pipeBase
from astropy import units as u
from astropy.table import Table
from lsst.daf.base import DateTime
import lsst.pex.config
from lsst.pipe.tasks.postprocess import TableVStack
from lsst.solsys.pipe.utils import df2numpy

_LOG = logging.getLogger(__name__)


class SubmitToMPCConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "day_obs"),
):
    consolidatedPureLinkages = pipeBase.connectionTypes.Input(
        doc="Input Source Tables to be concatenated",
        name="consolidatedPureLinkages",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "ssp_hypothesis_table", "day_obs"),
        multiple=True,
        minimum=0
    )
    discoveryVisitSummaries = pipeBase.connectionTypes.Input(
        doc="",
        name="visit_summary_dayobs_14",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "day_obs")
    )
    associatedSsSources = pipeBase.connectionTypes.Input(
        doc="",
        name="",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        multiple=True,
        minimum=0
    )
    associatedSsDiaSources = pipeBase.connectionTypes.Input( #Figure out how to not make this take 11 hours
        doc="",
        name="",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        multiple=True,
        deferLoad=True,
        minimum=0
    )
    # + Whatever we need to get seeing
    submittedAssociations = connectionTypes.Output(
        doc="Detections submitted to the MPC as associated asteroids",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_submitted_associations",
    )
    submittedDiscoveries = connectionTypes.Output(
        doc="Detections submitted to the MPC as new discoveries",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_submitted_discoveries",
    )
    submittedADESFiles = connectionTypes.Output(
        doc="The text of all ADES submissions to the MPC",
        dimensions=["instrument", "day_obs"],
        storageClass="ArrowAstropy",
        name="ssp_ADES_submissions",
    )


class SubmitToMPCConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=SubmitToMPCConnections):
    submissionDatabase = lsst.pex.config.Field(
        dtype=str,
        default='',
        doc='Location of the database containing our previous MPC submissions',
    )
    observatoryCode = lsst.pex.config.Field(
        dtype=str,
        default='X05',
        doc="Observatory MPC code, defaults to X05"
    )
    observatoryName = lsst.pex.config.Field(
        dtype=str,
        default='Vera C. Rubin Observatory',
        doc="Observatory name, defaults to Vera C. Rubin Observatory",
    )
    observerNames = lsst.pex.config.Field(
        dtype=str,
        default='',
        doc="Comma-separated list of observers (First initial, surname)",
    )
    measurerNames = lsst.pex.config.Field(
        dtype=str,
        default='',
        doc="Comma-separated list of measurers (First initial, surname)",
    )
    telescopeName = lsst.pex.config.Field(
        dtype=str,
        default='Simonyi Survey Telescope',
        doc="Name of telescope",
    )
    telescopeDesign = lsst.pex.config.Field(
        dtype=str,
        default='TODO: MJURIC',
        doc="Design of telescope",
    )
    telescopeAperture = lsst.pex.config.Field(
        dtype=float,
        default=-1,
        doc="Aperture of telescope in meters",
    )
    cameraDetector = lsst.pex.config.Field(
        dtype=str,
        default='CCD',
        doc='Type of detectors in camera',
    )
    softwareDiscovery = lsst.pex.config.Field(
        dtype=str,
        default='HelioLinC3D',
        doc='Software used for object discovery',
    )
    observatoryCode = lsst.pex.config.Field(
        dtype=str,
        default='X05',
        doc="Observatory code MPC, defaults to X05"
    )
    # TODO MJURIC: fRatio, filter, arraySize, pixelScale, coinvestigators, collaborators,
    #              fundingSource, comment, software for astrometry and fitOrder and photometry


class SubmitToMPCTask(pipeBase.PipelineTask):
    """Concatenate `sourceTable` list into a per-dayobs `sourceTable_dayobs`"""

    ConfigClass = SubmitToMPCConfig
    _DefaultName = "submitToMPC"

    def run(self, consolidatedPureLinkages, discoveryVisitSummaries, associatedSsSources,
            associatedSsDiaSources):
        """doc string
           here
        """
        ## Check MPC status
        # Read from MPC
        # Update our database to match MPC

        ## Submit new discoveries
        # Input consolidatedPureLinkages across many ssp_hypothesis_tables
        # De-duplicate them and check them against our previous MPC submission database
        # Prepare submission to MPC
        # Submit to MPC

        ## Submit new associations
        # Input ss_associated_sources across many (visit, detector)s
        # Check them against our previous MPC submission database
        # Prepare submission to MPC
        # Submit to MPC
        
        # Log metrics
        _LOG.info(f'Submitted X new objects and Y new associations to MPC. Refrained from submitting X overlapping.')
        return lsst.pipe.base.Struct(submittedAssociations=,
                                     submittedDiscoveries=,
                                     submittedADESFiles=
                                     )

