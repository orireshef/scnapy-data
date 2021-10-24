import os
import abc
from configobj import ConfigObj

class IScConfig():
    """
    An interface for scanpy steps configurations
    """

    @abc.abstractstaticmethod
    def get_seurat_defaults():
        pass

    @staticmethod
    def read_cfg_file(path: str, step: str) -> dict:
        try:
            assert os.path.isfile(path)
            config_obj = ConfigObj(path)
            return config_obj[step]
        except Exception as e:
            print(f"""
            Could not read config from file in {path}.
            Make sure the file and the step are defined correctly.
            """)

            
class PpConfig(IScConfig):
    """
    Preprocessing step configurations
    """

    def __init__(self, highest_expr_genes: dict,
    filter_cells: dict,
    filter_genes: dict,
    normalize_total: dict,
    highly_variable_genes: dict) -> None:
        self.highest_expr_genes = highest_expr_genes
        self.filter_cells = filter_cells
        self.filter_genes = filter_genes
        self.normalize_total = normalize_total
        self.highly_variable_genes = highly_variable_genes

    @staticmethod
    def get_seurat_defaults():
        return {
            "highest_expr_genes": {"n_top": 20},
            "filter_cells": {"min_genes": 200},
            "filter_genes": {"min_cells": 3},
            "normalize_total": {"target_sum": 1e4},
            "highly_variable_genes": {
                "min_mean": 0.0125,
                "max_mean": 3,
                "min_disp": 0.5
                }
        }