"""
Python reader to convert TrackmateXML to numpy
Version 1.0.0
(c) R.Harkes 2022 NKI

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import io
import json
import os
import logging
import sys
from contextlib import AbstractContextManager
from pathlib import Path
from types import TracebackType
from typing import Any, Union, Optional, Type, BinaryIO, List, Dict
import numpy as np
from lxml import etree
from dataclasses import dataclass


@dataclass
class AnalyzedTrack:
    parent: int
    cell: int
    spotids: np.ndarray[Any, np.dtype[np.int32]]
    track: bool


class TrackmateXMLFile(AbstractContextManager[Any]):
    """
    Context manager for handeling TrackmateXML files
    """

    def __init__(self, pth: Union[str, os.PathLike[Any]]) -> None:
        pth = Path(pth)
        if not pth.exists() or not pth.is_file():
            raise FileNotFoundError(pth)
        self.pth = pth

    def __enter__(self) -> "TrackmateXML":
        self.tmxml = TrackmateXML()
        self.file = io.FileIO(self.pth, "r")
        self.tmxml.loadstream(self.file)
        return self.tmxml

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        self.file.close()


class TrackmateXML:
    """
    Class to import a TrackmateXML. (c) R.Harkes GPLv3

    Please note that TrackMate is available through Fiji, and is based on a publication.
    If you use it successfully for your research please be so kind to cite the work:
    Tinevez, JY.; Perry, N. & Schindelin, J. et al. (2017), 'TrackMate: An open and extensible platform for single-particle tracking.', Methods 115: 80-90, PMID 27713081.
    https://www.ncbi.nlm.nih.gov/pubmed/27713081
    https://scholar.google.com/scholar?cluster=9846627681021220605
    """

    def __init__(self, loglevel: int = logging.ERROR) -> None:
        self.logger = logging.getLogger(__name__)
        if loglevel < logging.ERROR:
            handler = logging.StreamHandler(
                sys.stdout
            )  # set output to stdout instead of stderror
            handler.setLevel(loglevel)
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
        else:
            self.logger.setLevel(loglevel)
        self.logger.info("Logger Initialized")
        self.version = ""  # type: str
        self.spatialunits = ""  # type: str
        self.timeunits = ""  # type: str
        self.spotheader = []  # type:List[str]
        self.spots = np.zeros((0, 0), dtype=float)
        self.tracks = []  # type: List[np.ndarray[Any, np.dtype[np.int32]]]
        self.tracknames = []  # type: List[str]
        self.displaysettings = {}  # type: Dict[str, str]

    def __bool__(self) -> bool:
        if self.version:
            return True
        else:
            return False

    def loadstream(self, fp: BinaryIO) -> None:
        self._load(etree.parse(fp))

    def loadfile(self, pth: Union[str, os.PathLike[Any]]) -> None:
        """
        Load a TrackMate XML-file
        """
        self._load(etree.parse(pth))

    def loadtree(self, tree: etree._ElementTree) -> None:
        """
        Load a Trackmate XML-tree
        """
        self._load(tree)

    def gettraces(
        self,
        trackname: str,
        spot_property: str,
        duplicate_split: bool = False,
        break_split: bool = False,
    ) -> List[np.ndarray[Any, np.dtype[np.float64]]]:
        """
        Get traces from a trackname and a spot-property
        """
        tracks = self.analysetrack(trackname, duplicate_split, break_split)
        return [self.getproperty(track.spotids, spot_property) for track in tracks]

    def getproperty(
        self, spotids: np.ndarray[Any, np.dtype[np.int32]], spot_property: str
    ) -> np.ndarray[Any, np.dtype[np.float64]]:
        """
        Get properties from spotids
        """
        if spot_property not in self.spotheader:
            self.logger.error(f"{spot_property} not in spot properties")
            return np.zeros((0, 0), dtype=float)
        prop_idx = self.spotheader.index(spot_property)
        spotid_idx = self.spotheader.index("ID")
        res = np.zeros(len(spotids), dtype=np.float64)
        for i, s in enumerate(spotids):
            res[i] = self.spots[self.spots[:, spotid_idx] == s, prop_idx][0]
        return res

    def _load(self, tree: etree._ElementTree) -> None:
        self.logger.info("Getting version")
        self._getversion(tree.getroot())
        self.logger.info("Analyzing tree")
        self._analysetree(tree)
        self.logger.info("Finished loading")

    def getversion(self) -> str:
        """
        Get the version of the TrackmateXML data as string.
        """
        return self.version if self.version else ""

    def _getversion(self, root: etree._Element) -> None:
        if root.tag == "TrackMate":
            v = root.attrib.get("version", "")
            if v == "":
                self.logger.error(f"Invalid Version")
            elif type(v) == bytes:
                self.logger.error(f"Invalid Version")
            else:
                self.version = str(v)
        else:
            self.logger.error("Not a TrackMateXML")

    def _analysetree(self, tree: etree._ElementTree) -> None:
        root = tree.getroot()
        for element in root:
            if element.tag == "Log":
                self._getlog(element)
            elif element.tag == "Model":
                self._getmodel(element)
            elif element.tag == "Settings":
                self._getsettings(element)
            elif element.tag == "GUIState":
                self._get_gui_state(element)
            elif element.tag == "DisplaySettings":
                self._get_display_settings(element)
            else:
                self.logger.error(f"Unrecognised element {element}")

    def _getlog(self, element: etree._Element) -> None:
        self.log = element.text

    def _get_display_settings(self, element: etree._Element) -> None:
        self.displaysettings = json.loads(str(element.text))

    def _getmodel(self, element: etree._Element) -> None:
        self.spatialunits = str(element.attrib.get("spatialunits", ""))
        self.timeunits = str(element.attrib.get("timeunits", ""))
        for subelement in element:
            if subelement.tag == "FeatureDeclarations":
                pass  # would be nice if the feature declaration would actually list the features in the xml, but instead it lists all possible features
            elif subelement.tag == "AllSpots":
                self._getspots(subelement)
            elif subelement.tag == "AllTracks":
                self._gettracks(subelement)
            elif subelement.tag == "FilteredTracks":
                self._getfilteredtracks(subelement)
            else:
                self.logger.error(f"Unrecognised element {element}")

    def _getsettings(self, element: etree._Element) -> None:
        """
        Could be added, but is not required for displaying intensity tracks
        """
        pass

    def _getfilteredtracks(self, element: etree._Element) -> None:
        """
        Could be added, but is not required for displaying intensity tracks
        """
        pass

    def _get_gui_state(self, element: etree._Element) -> None:
        """
        Could be added, but does not seem usefull in python.
        """
        pass

    def _getspots(self, element: etree._Element) -> None:
        """
        Put all numeric spot data in a numpy array.
        """
        nspots = int(str(element.attrib.get("nspots", "0")))
        # construct header
        spotid = 0
        spot = element[0][0]
        keys = [str(a) for a in spot.attrib]
        for k in keys:
            try:
                float(spot.attrib[k])
            except ValueError:  # remove keys we cannot convert to floats
                keys.remove(k)
        self.spotheader = keys
        self.spots = np.zeros((nspots, len(keys)), dtype=np.float64)
        for sif in element:
            for spot in sif:
                for i, k in enumerate(keys):
                    self.spots[spotid, i] = float(str(spot.attrib.get(k, "nan")))
                spotid += 1

    def _gettracks(self, element: etree._Element) -> None:
        """
        Importing only source and target into a numpy array and listing the trackname.
        We could import more, but it is not needed for displaying intensity tracks.
        """
        for track in element:
            t = np.zeros((len(track), 2), dtype=np.int32)
            for i, edge in enumerate(track):
                t[i, 0] = int(str(edge.attrib.get("SPOT_SOURCE_ID", "-1")))
                t[i, 1] = int(str(edge.attrib.get("SPOT_TARGET_ID", "-1")))
            self.tracks.append(t)
            self.tracknames.append(str(track.attrib.get("name", "unknown")))

    def analysetrack(
        self, trackname: str, duplicate_split: bool = False, break_split: bool = False
    ) -> List[AnalyzedTrack]:
        """
        Traces a track to find the sequence of spotids
        """
        trackid = self.tracknames.index(trackname)
        return self.analysetrackid(trackid, duplicate_split, break_split)

    def analysetrackid(
        self, trackid: int, duplicate_split: bool = False, break_split: bool = False
    ) -> List[AnalyzedTrack]:
        """
        Traces a track to find the sequence of spotids
        """
        track = self.tracks[trackid]
        unique_sources = np.setdiff1d(track[:, 0], track[:, 1], assume_unique=True)
        if unique_sources.size != 1:
            self.logger.error(
                f"Track has {unique_sources.size} startingpoints. Cannot follow track."
            )
            return []
        cellid = 1
        traced_tracks = [
            AnalyzedTrack(
                parent=0,
                cell=cellid,
                spotids=track[np.argwhere(track[:, 0] == unique_sources), :].flatten(),
                track=True,
            )
        ]
        while any([x.track for x in traced_tracks]):
            for traced_track in traced_tracks:
                if not traced_track.track:
                    continue
                traced_track_ids = traced_track.spotids
                targetidx = np.argwhere(
                    track[:, 0] == traced_track_ids[-1]
                ).flatten()  # do we find the last index in the sources?
                if targetidx.size == 0:
                    traced_track.track = False  # reached the end of the track
                elif targetidx.size == 1:
                    traced_track.spotids = np.concatenate(
                        (traced_track_ids, track[targetidx, 1])
                    )  # append target to track
                else:  # multiple targets, a split
                    if duplicate_split:  # a copy of the history is added to each track
                        for i in range(1, targetidx.size):
                            spotids = np.concatenate(
                                (traced_track_ids, track[targetidx[i], 1].flatten())
                            )
                            cellid += 1
                            new_track = AnalyzedTrack(
                                parent=traced_track.cell,
                                cell=cellid,
                                spotids=spotids,
                                track=True,
                            )
                            traced_tracks.append(new_track)
                    else:
                        for i in range(1, targetidx.size):
                            spotids = track[targetidx[i], :]
                            cellid += 1
                            new_track = AnalyzedTrack(
                                parent=traced_track.cell,
                                cell=cellid,
                                spotids=spotids,
                                track=True,
                            )
                            traced_tracks.append(new_track)
                    if break_split:
                        traced_track.track = False  # end the parent
                        if duplicate_split:
                            spotids = np.concatenate(
                                (traced_track_ids, track[targetidx[0], 1].flatten())
                            )
                        else:
                            spotids = track[targetidx[0], :]  # start child
                        cellid += 1
                        new_track = AnalyzedTrack(
                            parent=traced_track.cell,
                            cell=cellid,
                            spotids=spotids,
                            track=True,
                        )
                        traced_tracks.append(new_track)
                    else:
                        traced_track.spotids = np.concatenate(
                            (traced_track_ids, track[targetidx[0], 1].flatten())
                        )  # append target to track
        return traced_tracks