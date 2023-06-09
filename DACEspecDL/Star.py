from __future__ import annotations

from astroquery.simbad import Simbad

import shutil
from dataclasses import dataclass
from functools import cached_property, lru_cache
from pathlib import Path
from typing import List, Iterable, Any, Dict, Optional

from dace_query import Dace
from dace_query.spectroscopy import Spectroscopy
import matplotlib.pyplot as plt
from loguru import logger


@dataclass
class Star:
    name: str
    pipe_KW: Optional[Dict[str, str]] = None
    api_user: Optional[str] = None

    @cached_property
    def _data(self):
        """
        subInstrument
            Pipeline level
                OBS mode
        :return:
        """
        if self.api_user is not None:
            desired_key = Dace._DaceClass__dace_rc_config.items('amiguel')[0][1]
            Dace._DaceClass__dace_rc_config.set('user', 'key', desired_key)

        return Spectroscopy.get_timeseries(target=self.name,
                                           sorted_by_instrument=True,
                                           )

    def get_available_instruments(self) -> List[str]:
        return list(self._data.keys())

    def get_pipelines_of_instrument(self, instrument: str) -> List[str]:
        return list(self._data[instrument].keys())

    def get_OBS_modes(self, instrument: str, pipeline: str) -> List[str]:
        """
        Get the modes of observations
        :param instrument:
        :param pipeline:
        :return:
        """
        if instrument not in self.get_available_instruments() or pipeline not in self.get_pipelines_of_instrument(instrument):
            raise Exception(f"Either the instrument or the pipeline don't exist - {instrument}::{pipeline}")
        return list(self._data[instrument][pipeline].keys())

    def get_header_info(self,
                        kw_to_get: List[str],
                        instrument: Optional[str] = None,
                        OBS_mode: Optional[str] = None,
                        pipe_identifier: Optional[str] = None,
                        ) -> Dict[str, Any]:

        data_dict = {}
        for (instrument_name, pipe_name, OBS_name, data) in self.data_to_iterate_over(metric_name="rv",
                                                                                      pipe_identifier=pipe_identifier,
                                                                                      OBS_mode=OBS_mode,
                                                                                      instrument=instrument,
                                                                                      get_full_dict=True
                                                                                      ):
            if instrument_name not in data_dict:
                data_dict[instrument_name] = {}

            if OBS_name not in data_dict[instrument_name]:
                data_dict[instrument_name] = {OBS_name: {}}

            if pipe_name not in data_dict[instrument_name][OBS_name]:
                data_dict[instrument_name][OBS_name][pipe_name] = {i: data[i] for i in kw_to_get}

        return data_dict

    def get_RVs(self, instrument: Optional[str] = None,
                OBS_mode: Optional[str] = None,
                pipe_identifier: Optional[str] = None
                ) -> Dict[str, float]:
        return self.get_header_info(kw_to_get=["rjd", "rv", "rv_err"],
                                    instrument=instrument,
                                    OBS_mode=OBS_mode,
                                    pipe_identifier=pipe_identifier
                                    )

    def download_data(self,
                      output_path: Path,
                      force_download: bool = False,
                      instrument: Optional[str] = None,
                      OBS_mode: Optional[str] = None,
                      pipe_identifier: Optional[str] = None,
                      file_type: str = "all",
                      common_root_folder: bool = True,  # for the unpack
                      unzip: bool = True,  #
                      allow_subfolders: bool = True
                      ) -> None:
        """

        :param output_path: Path in which the spectra will be stored
        :param force_download: Download, even if the files exist on disk
        :param instrument: sequence of characters that exist in the instrument name
        :param OBS_mode: sequence of characters that exist in the OBS mode
        :param pipe_identifier: sequence of characters that exist in the pipeline name
        :param file_type:
        :param common_root_folder: After unzipping the data, keep everything on the created subfolder
        :param unzip: Trigger unzip of the downloaded .tar file
        :param allow_subfolders: Allow creation of subfolders to divide by instrument and pipeline name
        :return:
        """
        count = 0
        file_list_to_download = {}

        tar_name = "result.tar.gz"
        if instrument.upper() == "HARPSN":
            logger.warning("The name of HARPSN on DACE is HARPN; Updating the name")
            instrument = "HARPN"

        for (*_, item) in self.data_to_iterate_over(metric_name="raw_file",
                                                    instrument=instrument, OBS_mode=OBS_mode,
                                                    pipe_identifier=pipe_identifier
                                                    ):
            count += 1
            for file in item:
                store_inst_name, store_pipe_name, *_ = file.split("/")

                if allow_subfolders:
                    modifier = output_path / store_inst_name / store_pipe_name
                else:
                    modifier = output_path

                if modifier not in file_list_to_download:
                    file_list_to_download[modifier] = []

                similar_files = list(modifier.glob(f'**/{file.split("/")[-1].split(".fits")[0]}*'))
                if len(similar_files) != 0 and not force_download:
                    # The file already exists on disk
                    continue
                file_list_to_download[modifier].append(file)

        logger.info("Launching downloads")
        for disk_path, filelist in file_list_to_download.items():
            if len(filelist) == 0:
                logger.warning("All files already exist in the specific disk location")
                continue
            logger.info(f"Triggering the download of {len(filelist)} files")

            disk_path.mkdir(exist_ok=True, parents=True)
            Spectroscopy.download_files(files=filelist,
                                        output_directory=disk_path,
                                        file_type=file_type,
                                        output_filename=tar_name
                                        )

            if unzip:
                import tarfile
                with tarfile.open(disk_path / tar_name) as tar:
                    tar.extractall(path=disk_path)

                if common_root_folder:
                    all_paths = []
                    for i in disk_path.iterdir():
                        if not i.is_dir():
                            continue
                        all_paths.extend(list(i.glob("*.fits")))

                    store_path = disk_path.as_posix()

                    for path in all_paths:
                        shutil.move(path.as_posix(), store_path)

                    for path in disk_path.iterdir():
                        if not path.is_dir():
                            continue

                        path.rmdir()

        print("Total sum: ", count)

    def data_to_iterate_over(self,
                             metric_name: str,
                             instrument: Optional[str] = None,
                             OBS_mode: Optional[str] = None,
                             pipe_identifier: Optional[str] = None,
                             get_full_dict=False
                             ):
        """

        :param metric_name:
        :param instrument:
        :param OBS_mode:
        :param pipe_identifier:
        """
        reject_iter = lambda x, y: x is not None and x not in y

        for instrument_name, data_1 in self._data.items():
            if reject_iter(instrument, instrument_name):
                continue

            for pipe_name, data_2 in data_1.items():
                if reject_iter(pipe_identifier, pipe_name):
                    continue

                for OBS_name, data_3 in data_2.items():
                    if reject_iter(OBS_mode, OBS_name):
                        continue
                    if get_full_dict:
                        yield instrument_name, pipe_name, OBS_name, data_3
                    else:
                        yield instrument_name, pipe_name, OBS_name, data_3[metric_name]

    def get_metrics_of_instrument(self, instrument: str, metric: str | Iterable[str]) -> Dict[str, List[Any]]:
        """
        TODO: add masking of OBS based on the mode
        TODO: allow regex for the instrument name

        :param instrument:
        :param metric:
        :return:
        """

        if not isinstance(metric, Iterable):
            metric = [metric]

        pipe_KW = self.get_pipe_KW(instrument)
        data: Dict[str, List] = {i: [] for i in metric}

        for pipe_name, data in self._data[instrument]:
            if pipe_KW not in pipe_name:
                # Do this to ensure that this doesn't break in ESPRESSO, due to the HR names also
                # being in the pipe name
                continue
            for OBS_name, data in data.items():  # Iterate over all possible OBS modes
                for met in metric:
                    data[met].append(data[met])
        return data

    def plot_Rvs(self):
        fig, axis = plt.subplots()

    def __str__(self):
        data = self._data

        for key, value in data.items():
            print("High level key: ", key)
            for key, value in value.items():
                print("Mid level key:", key)
                for key, value in value.items():
                    print("low level key:", key)

        return self.name

    def get_pipe_KW(self, instrument):
        available_pipes = self.get_pipelines_of_instrument(instrument)

        if len(available_pipes) == 1:
            return available_pipes[0]

        try:
            if self.pipe_KW is not None:
                return self.pipe_KW[instrument]
        except KeyError:
            logger.warning(f"User provided the desired pipe-KW dictionary, but it does not specify pipe for {instrument}")

        if "ESPRESSO" in instrument:
            return "3.0.0"
        elif "HARPN" in instrument:
            return "2.3.5"

    @property
    def spectral_type(self) -> str:
        """
        Get the spectral type from SIMBAD
        :return:
        """

        aSimbad = Simbad()
        aSimbad.add_votable_fields("sptype")
        r = aSimbad.query_object(self.name)

        print(r.keys())
        if r is None:
            raise Exception(f"Couldn't find the object {self.name} in SIMBAD")
        ST = r["SP_TYPE"][0][:2].replace("d", "")
        return ST

    @property
    def aliases(self):
        return Simbad.query_objectids(self.name)
