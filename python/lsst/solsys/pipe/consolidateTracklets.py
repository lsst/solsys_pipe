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

"""A task which takes per-dayobs source tables and visit summaries and
consolidates them into multi-day tables. 
"""

__all__ = [
    "ConsolidateTrackletsTask",
    "ConsolidateTrackletsConfig",
]


import logging
import warnings
import lsst.pex.config
import lsst.pipe.base as pipeBase
import numpy as np
from . import utils
from astropy import units as u
import astropy.table as tb
from heliolinx import solarsyst_dyn_geo as solardg
from lsst.daf.base import DateTime
from lsst.pipe.tasks.postprocess import TableVStack
from lsst.resources import ResourcePath
from lsst.solsys.pipe.utils import df2numpy

_LOG = logging.getLogger(__name__)
warnings.filterwarnings("ignore")


class ConsolidateTrackletsConnections(
    pipeBase.PipelineTaskConnections,
    defaultTemplates={
        "diaSourceInputName": "dia_source_dayobs",
        "earthStateInputName": "sspEarthState",
        "diaSourceOutputName": "dia_source_dayobs_14",
        "visitSummaryInputName": "visit_summary_dayobs",
        "visitInfoOutputName": "visit_summary_dayobs_14",
    },
    dimensions=("instrument", "day_obs", "ssp_hypothesis_table"),
):
    inputVisitSummaries = pipeBase.connectionTypes.Input(
        doc="Per-dayobs consolidated exposure metadata",
        name="{visitSummaryInputName}",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "day_obs"),
        multiple=True,
    )
    inputTrackletSources = pipeBase.connectionTypes.Input(
        doc="sources that got included in tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_source_dayobs",
        multiple=True,
    )
    inputTracklets = pipeBase.connectionTypes.Input(
        doc="summary data for tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_dayobs",
        multiple=True,
    )
    inputTrackletToSource = pipeBase.connectionTypes.Input(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_to_source_dayobs",
        multiple=True,
    )
    outputTrackletSources = pipeBase.connectionTypes.Output(
        doc="sources that got included in tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_source_dayobs_14"
    )
    outputTracklets = pipeBase.connectionTypes.Output(
        doc="summary data for tracklets",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_dayobs_14",
    )
    outputTrackletToSource = pipeBase.connectionTypes.Output(
        doc="indices connecting tracklets to sspTrackletSources",
        dimensions=["instrument", "day_obs", "ssp_hypothesis_table"],
        storageClass="ArrowAstropy",
        name="ssp_tracklet_to_source_dayobs_14"
    )
    outputVisitSummaries = pipeBase.connectionTypes.Output(
        doc="Concatenated visit summary from all day_obs in the input with day_obs"
        "dimension of the latest day_obs in the input.",
        name="visit_summary_dayobs_14",
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
            inputs = adjuster.get_inputs(data_id)
            # In the three for loops below, we loop over all input refs with
            # the same data_id as the current data_id in the loop and add
            # their input data to the quantum for the latest day_obs.
            input_data_types = ["inputVisitSummaries", "inputTracklets",
                                "inputTrackletToSource", "inputTrackletSources"]
            for input_data_type in input_data_types:
                for input_data_id in inputs[input_data_type]:
                    adjuster.add_input(data_id_latest, input_data_type, input_data_id)
            adjuster.remove_quantum(data_id)

        # Log that the last day_obs is being kept as a reference.
        _LOG.info(
            f"Combined inputs from {len(to_do)} quanta into one quantum "
            f"under reference day_obs {data_id_latest['day_obs']}."
        )


class ConsolidateTrackletsConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=ConsolidateTrackletsConnections
):
    pass


class ConsolidateTrackletsTask(pipeBase.PipelineTask):
    """Concatenate `sourceTable` list into a per-dayobs `sourceTable_dayobs`"""

    ConfigClass = ConsolidateTrackletsConfig
    _DefaultName = "consolidateTracklets"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Concatenate the input source tables into a single output table."""

        # Let's make the order deterministic by sorting the input Refs.
        inputRefs.inputVisitSummaries.sort(key=lambda x: (x.dataId["day_obs"]))
        inputRefs.inputTrackletSources.sort(key=lambda x: (x.dataId["day_obs"]))
        inputRefs.inputTracklets.sort(key=lambda x: (x.dataId["day_obs"]))
        inputRefs.inputTrackletToSource.sort(key=lambda x: (x.dataId["day_obs"]))

        inputs = butlerQC.get(inputRefs)
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputVisitSummaries, inputTrackletSources, inputTracklets, inputTrackletToSource):
        n_tables = [len(inputVisitSummaries), len(inputTrackletSources),
                        len(inputTracklets), len(inputTrackletToSource)]
        assert min(n_tables) == max(n_tables), "Not all lists of tables same length!"
        n_tables = n_tables[0]
        if n_tables < 1:
            raise pipeBase.NoWorkFound
        n_tracklets, n_sources, n_visits = 0, 0, 0  # cumulative numbers of 
        for i in range(n_tables):
            # no need to update inputVisitSummaries[i]
            # update inputTrackletSources[i]
            print('\n\n\n')
            print(max(inputTrackletToSource[i]['i1']) + 1, len(inputTracklets[i]))
            print(max(inputTrackletToSource[i]['i2']) + 1, len(inputTrackletSources[i]))
            assert max(inputTrackletToSource[i]['i1']) + 1 == len(inputTracklets[i]), "shuffled tables?!"
            assert max(inputTrackletToSource[i]['i2']) + 1 == len(inputTrackletSources[i]), "shuffled tables?!"
            inputTrackletSources[i]['index'] += n_sources
            inputTrackletSources[i]['image'] += n_visits
            # update inputTracklets[i]
            inputTracklets[i]['Img1'] += n_visits
            inputTracklets[i]['Img2'] += n_visits
            inputTracklets[i]['trk_ID'] += n_tracklets
            # update inputTrackletToSource[i]
            inputTrackletToSource[i]['i1'] += n_tracklets
            inputTrackletToSource[i]['i2'] += n_sources
            # update cumulative 
            n_tracklets += len(inputTracklets[i])
            n_sources += len(inputTrackletSources[i])
            n_visits += len(inputVisitSummaries[i])
        visitSummary = tb.vstack(inputVisitSummaries)
        trackletSources = tb.vstack(inputTrackletSources)
        tracklets = tb.vstack(inputTracklets)
        trackletToSource = tb.vstack(inputTrackletToSource)

        return lsst.pipe.base.Struct(
            outputVisitSummaries=visitSummary, outputTrackletSources=trackletSources, outputTracklets=tracklets,
            outputTrackletToSource=trackletToSource)
