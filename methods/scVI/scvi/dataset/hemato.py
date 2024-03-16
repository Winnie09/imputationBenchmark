import os
import logging
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class HematoDataset(DownloadableDataset):
    """Loads the hemato dataset.

    This dataset contains continuous gene expression variations from hematopoeitic progenitor cells [31] contains
    4,016 cells and 7,397 genes. We removed the library basal-bm1 which was of poor quality based on authors
    recommendation. We use their population balance analysis result as a potential function for differentiation.

    Examples:
        >>> gene_dataset = HematoDataset()
    """

    def __init__(
        self, save_path: str = "data/HEMATO/", delayed_populating: bool = False
    ):
        self.gene_names_filename = "bBM.filtered_gene_list.paper.txt"
        self.spring_and_pba_filename = "bBM.spring_and_pba.csv"
        self.cell_types_levels = [
            "Erythroid",
            "Granulocytic Neutrophil",
            "Lymphocytic",
            "Dendritic",
            "Megakaryocytic",
            "Monocytic",
            "Basophilic",
        ]
        super().__init__(
            urls=[
                "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2388072&format=file&"
                "file=GSM2388072%5Fbasal%5Fbone%5Fmarrow%2Eraw%5Fumifm%5Fcounts%2Ecsv%2Egz",
                "https://github.com/romain-lopez/scVI-reproducibility/raw/master/additional/data.zip",
            ],
            filenames=["bBM.raw_umifm_counts.csv.gz", "data.zip"],
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing Hemato data")

        if len(os.listdir(self.save_path)) == 2:  # nothing extracted yet
            with ZipFile(os.path.join(self.save_path, "data.zip"), "r") as zip:
                zip.extractall(path=Path(self.save_path).parent)
        raw_counts = pd.read_csv(
            os.path.join(self.save_path, self.filenames[0]), compression="gzip"
        )

        # remove this library to avoid dealing with batch effects
        raw_counts.drop(
            raw_counts.index[raw_counts["library_id"] == "basal_bm1"], inplace=True
        )

        spring_and_pba = pd.read_csv(
            os.path.join(self.save_path, self.spring_and_pba_filename)
        )
        gene_names = np.loadtxt(
            os.path.join(self.save_path, self.gene_names_filename), dtype=np.str
        )

        data = raw_counts.merge(spring_and_pba, how="inner")
        expression_data = data[gene_names]
        x_spring = data["x_spring"].values
        y_spring = data["y_spring"].values

        self.meta = data[
            ["Potential", "Pr_Er", "Pr_Gr", "Pr_Ly", "Pr_DC", "Pr_Mk", "Pr_Mo", "Pr_Ba"]
        ]

        def logit(p):
            p = np.copy(p.values)
            p[p == 0] = np.min(p[p > 0])
            p[p == 1] = np.max(p[p < 1])
            return np.log(p / (1 - p))

        labels = logit(self.meta.iloc[:, 2]) - logit(self.meta.iloc[:, 1])
        expression_data = expression_data.values

        logger.info("Finished preprocessing Hemato data")
        self.populate_from_data(
            X=expression_data,
            labels=labels,
            gene_names=gene_names,
            cell_attributes_dict={"x_coord": x_spring, "y_coord": y_spring},
        )
        self.filter_cells_by_count()
