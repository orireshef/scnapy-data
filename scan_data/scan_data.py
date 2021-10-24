import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from .sc_config import PpConfig


class ScanData(ad):
    """
    Extends anndata object with scanpy functionality
    """

    sc_config: PpConfig

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
        self.read_data(path, **kwargs)
        