import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from .sc_step import PpStep
from configobj import ConfigObj


class ScanData():
    """
    Extends anndata object with scanpy functionality
    """

    _sc_config: ConfigObj
    _is_configured = False

    @staticmethod
    def read_cfg_file(path: str) -> ConfigObj:
        try:
            assert os.path.isfile(path)
            config_obj = ConfigObj(path)
            return config_obj
        except Exception as e:
            print(f"""
            Could not read config from file in {path}.
            Make sure the file and the step are defined correctly.
            """)
            raise(e)

    # define a method to create a default sc_config obj (seurat)
    
    def read_data(self, path: str, **kwargs):
        try:
            if os.path.isdir:
                print(f"attempting to read from dir: {path}")
                self.adata = sc.read_10x_mtx(path, **kwargs)
            else:
                print(f"attempting to read from file: {path}")
                self.adata = sc.read(path, **kwargs)

        except Exception as e:
            print(f"Failed to load {path} into an anndata object.")
            raise e

    def __init__(self, path: str, sc_config = None, **kwargs):
        """
        Instaciate a ScanData object from a file or directory
        containing a Data file.
        :param path: A path to a file or directory containing a data file
        :param sc_config: a path to a scanpy.cfg file or a ScConfig object. if not provided
        will use default seurat configurations.
        :param kwargs: parameters passed to read or read_10x_mtx
        """
        if isinstance(sc_config, str):
            sc_config = self.read_cfg_file(sc_config)
        elif sc_config is None:

        self.sc_config = sc_config
        self.read_data(path, **kwargs)

    @property
    def adata(self):
        return self._adata

    @adata.setter
    def adata(self, adata):
        self._adata = adata
    
    @property
    def sc_config(self):
        return self._sc_config

    @sc_config.setter
    def sc_config(self, sc_config: ConfigObj):
        if self._is_configured:
            raise Exception("Cannot reset sc config after ScanData created") 