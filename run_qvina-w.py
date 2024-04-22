import argparse
import os 
import subprocess
import warnings
import json

from meeko import MoleculePreparation
from moleculekit.molecule import Molecule 
from moleculekit.projections.metricsecondarystructure import MetricSecondaryStructure
import pandas as pd 
import numpy as np
import yaml
from rdkit import Chem 
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Geometry import Point3D
from tqdm import tqdm

from src.converters import *
from src.dockings import *
from src.fetchers import *
from src.loaders import *

def get_parser():
    parser = argparse.ArgumentParser(
        description="",
        usage=f"python {os.path.basename(__file__)} -c CONFIG_FILE"
    )
    parser.add_argument(
        "-c", "--config", type=str, required=True,
        help="path to a config file"
    )
    return parser.parse_args()

def set_default_config(conf):
    conf.setdefault("dataset", None)
    conf.setdefault("ligand_type", "sdf")
    conf.setdefault("protein_type", "pdb")
    conf.setdefault("docking_type", "specified")
    conf.setdefault("afdb_path", None)
    conf.setdefault("ligand_list", "ligand_list.csv")
    conf.setdefault("docking_pair", "docking_pair_list.csv")
    conf.setdefault("num_modes", 16)
    conf.setdefault("exhaustiveness", 16)
    conf.setdefault("output_dir", "../output")


def main():
    args = get_parser()
    with open(args.config, "r") as f:
        conf = yaml.load(f, Loader=yaml.SafeLoader)
    set_default_config(conf)
    
    # ディレクトリの作成
    os.makedirs("../input/ligand_dir", exist_ok=True)
    os.makedirs("../input/receptor_dir", exist_ok=True)
    os.makedirs(conf["output_dir"], exist_ok=True)

    print(f"========== Configuration ==========")
    for k, v in conf.items():
        print(f"{k}: {v}")
    print(f"===================================")
    
    print("Input Loading...")
    if conf['ligand_type'] == 'smiles':
        if conf["protein_type"] == "sequence":
            #  csvファイルを読み込む.
            csv_loader = CSVLoader("../input/%s"%conf["activity_file"])
            csv_loader.load()
        elif conf["protein_type"] == "pdb" or "pdbqt":
            # 何もしなくていい
            pass
        else:
            # エラー
            raise ValueError("protein_type must be sequence, pdb, or pdbqt.")
    elif conf['ligand_type'] == 'sdf':
        # リガンドリストの読み込み
        if conf["docking_type"] == "pairwise":
            list_loader = LigListLoader("../input/%s"%conf["ligand_list"])
            list_loader.load()
        elif conf["docking_type"] == "specified":
            list_loader = LigListLoader("../input/%s"%conf["ligand_list"])
            list_loader.load()
            docking_pair_loader = DockingPairLoader("../input/%s"%conf["docking_pair"])
            docking_pair_loader.load()
    elif conf["ligand_type"] == "pdbqt":
        pass
    else:
        raise ValueError("ligand_type must be smiles, sdf, or pdbqt.")

    print("Ligand Preparation...")
    # ligandの変換
    os.makedirs("../input/ligand_pdbqt", exist_ok=True)
    
    if conf['ligand_type'] == 'smiles':
        smiles_converter = SMILESConverter(csv_loader)
        smiles_converter.convert()
        
    elif conf['ligand_type'] == 'sdf':
        sdf_converter = SDFConverter(list_loader)
        sdf_converter.convert()
        
    elif conf["ligand_type"] == "pdbqt":
        pass
    
    else:
        raise ValueError("ligand_type must be smiles, sdf, or pdbqt.")
    
        
    print("Protein Preparation...")
    # proteinの変換
    os.makedirs("../input/receptor_pdbqt", exist_ok=True)
    if conf['protein_type'] == 'sequence':
        # AlphaFoldDBから構造を取得
        afdb_path = conf["afdb_path"]
        afdb_fetcher = AFDBFetcher(csv_loader, afdb_path)
        afdb_fetcher.fetch()
        # PDB -> PDBQT
        pdb_converter = PDBConverter(path="../input/receptor_dir", conf=conf)
        pdb_converter.convert()
        
    elif conf['protein_type'] == 'pdb':
        # PDB -> PDBQT
        if conf["docking_type"] == "specified":
            pdb_converter = PDBConverter(path="../input/receptor_dir", conf=conf, loader=docking_pair_loader)
        else:
            pdb_converter = PDBConverter(path="../input/receptor_dir", conf=conf)
        pdb_converter.convert()
    
    elif conf["protein_type"] == "pdbqt":
        pass
    
    else:
        raise ValueError("protein_type must be sequence, pdb, or pdbqt.")

"""
    print("Docking...")    
    if conf["docking_type"] == "pairwise":
        # pairwise docking
        docking = PairwiseDocking(conf, list_loader)
        docking.dock()
    
    elif conf["docking_type"] == "specified":
        # docking in specified pairs
        docking = SpecifiedDocking(conf, docking_pair_loader)
        docking.dock()
"""

if __name__ == "__main__":
    main()
