import os
import warnings
import json 

import numpy as np
import pandas as pd
from meeko import MoleculePreparation
from moleculekit.molecule import Molecule 
from moleculekit.projections.metricsecondarystructure import MetricSecondaryStructure
from rdkit import Chem 
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Geometry import Point3D
from tqdm import tqdm
import glob


class CSVLoader:
    def __init__(self, path):
        self.path = path

    def load(self):
        df_data = pd.read_csv(self.path)
        self.df = df_data
        self.uniprot_ids = pd.factorize(df_data["target_accession"])[1]
        if df_data["drug_id"].nunique() != df_data["smiles"].nunique():
            warnings.warn("drug_id is not appropriate. Please renumber the unique ids.", UserWarning)
        unique_drugs = df_data[["drug_id", "smiles"]].drop_duplicates(keep="first")
        self.drug_ids = unique_drugs["drug_id"].tolist()
        self.smiles = unique_drugs["smiles"].tolist()
        
class LigListLoader:
    def __init__(self, path):
        self.path = path

    def load(self):
        df_data = pd.read_csv(self.path)
        unique_drugs = df_data["ID"].unique()
        self.uids = list(unique_drugs)

class DockingPairLoader:
    def __init__(self, path):
        self.path = path
    
    def load(self):
        df_data = pd.read_csv(self.path)
        self.df = df_data
        self.ligand_name = df_data["ligand_name"].tolist()
        self.receptor_name = df_data["receptor_name"].tolist()