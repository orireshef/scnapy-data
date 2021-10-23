import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os


class ScanData():
    """
    Extends anndata object with scanpy functionality
    """

    def __init__(self, path: str, **kwargs):
        """
        Instaciate a ScanData object from a file or directory
        containing a Data file.
        :param path: A path to a file or directory containing a data file
        :param kwargs: parameters passed to read or read_10x_mtx
        """

        try:
            if os.path.isdir:
                print(f"attempting to read from dir: {path}")
                self.adata = sc.read_10x_mtx(path, kwargs)
            else:
                print(f"attempting to read from file: {path}")
                self.adata = sc.read(path, kwargs)
        except Exception as e:
            print(f"Failed to load {path} into an anndata object.")
            raise e