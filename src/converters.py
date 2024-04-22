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

class SMILESConverter:
    def __init__(self, loader):
        self.loader = loader
        
    def _convert(self, drug_id, smi):
        lig_name = f"drug_id-{drug_id}"
        pdbqt_lig = "../input/ligand_pdbqt/%s.pdbqt" % lig_name
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        #mol_conf = mol.GetConformer(-1)
        AllChem.UFFOptimizeMolecule(mol)
        #centroid = list(rdMolTransforms.ComputeCentroid(mol_conf))
        #vina_center = [0, 0, 0]
        #tr = [vina_center[i] - centroid[i] for i in range(3)]
        #for i, p in enumerate(mol_conf.GetPositions()):
        #    mol_conf.SetAtomPosition(i, Point3D(p[0]+tr[0], p[1]+tr[1], p[2]+tr[2]))
        #mol_prep = MoleculePreparation()
        writer = Chem.SDWriter("../input/ligand_dir/%s.sdf" % lig_name)
        mol.SetProp('smiles', smi)
        writer.write(mol)
        writer.close()
        #mol_prep.prepare(mol)
        #with open(pdbqt_lig, "w") as f:
        #    f.write(mol_prep.write_pdbqt_string())
        obabel = '/home/apps/openbabel/2.4.1/bin/obabel'
        sdf_lig_name = "../input/ligand_dir/%s.sdf" % lig_name
        pdbqt_lig_name = "../input/ligand_pdbqt/%s.pdbqt" % lig_name
        cmd = "%s -h -isdf %s -opdb -O %s" % (obabel, sdf_lig_name, pdbqt_lig_name)
        os.system(cmd)
        
    def convert(self):
        for drug_id, smi in tqdm(zip(self.loader.drug_ids, self.loader.smiles)):
            lig_name = f"drug_id-{drug_id}"
            pdbqt_lig_name = "../input/ligand_pdbqt/%s.pdbqt" % lig_name
            if not os.path.exists(pdbqt_lig_name):
                try:
                    self._convert(drug_id, smi)
                except Exception as e:
                    print(e)
                    print("cannot convert into pdbqt: drug_id = %s, smiles = %s" %(drug_id, smi))
                    pass

class PDBConverter:
    def __init__(self,path, conf, loader=None):
        self.path = path
        self.conf = conf
        self.loader = loader
    
    def _convert(self, rec_name):
        pdb_rec_path = os.path.join(self.path, "%s.pdb"%rec_name)
        mol = Molecule(pdb_rec_path)
        mol_cp = mol.copy()
        mol_cp.filter('protein')
        if self.conf["dataset"] == "pdbbind":
            mol_cp = self._select_chain(rec_name, mol_cp) 
        metr = MetricSecondaryStructure(simplified=False)
        labels = metr.project(mol_cp)
        mol_without_loop = mol_cp.copy()
        resid = " ".join(str(i) for i in np.where(labels==10)[1])
        mol_without_loop.remove("resid %s"%resid)
        mol_without_loop.center([0, 0, 0])
        coords = mol_without_loop.coords
        center = (coords.max(axis=0) + coords.min(axis=0))/2
        center = np.round(center, decimals=3)
        box_size = (coords.max(axis=0) -  coords.min(axis=0)) + 8
        box_size = np.round(box_size, decimals=3)
        
        with open(os.path.join("../input/receptor_pdbqt", "%s_box.json"%rec_name), "w") as f:
            bconf = {}
            bconf["center_x"] = str(center[0][0])
            bconf["center_y"] = str(center[1][0])
            bconf["center_z"] = str(center[2][0])
            bconf["size_x"] = str(box_size[0][0])
            bconf["size_y"] = str(box_size[1][0])
            bconf["size_z"] = str(box_size[2][0])
            json.dump(bconf, f, indent=4)
            
        pdb_rec_centered = os.path.join("../input/receptor_centered", "%s_centered.pdb"%rec_name)
        mol_cp.write(pdb_rec_centered)
        
        pdbqt_rec_path = "../input/receptor_pdbqt/%s.pdbqt"%rec_name
        adtools = '/home/mabiao/apps/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'
        cmd = "%s -r %s -o %s -A checkhydrogens" % (adtools, pdb_rec_centered, pdbqt_rec_path)
        os.system(cmd)
    
    def _select_chain(self, rec_name, mol):
        mol_cp = mol.copy()
        rec_name = rec_name.split("_protein")[0]
        pdb_pocket_path = os.path.join("../input/receptor_pocket", "%s_pocket.pdb"%(rec_name))
        mol_pocket = Molecule(pdb_pocket_path)
        mol_pocket.filter('protein')
        uniq_ch = np.unique(mol_pocket.get('chain'))
        
        maj_ch = uniq_ch[0]
        len_maj_ch = len(np.unique(mol_pocket.get('resid', 'chain %s'%maj_ch)))
        for ch in uniq_ch:
            uniq_res = np.unique(mol_pocket.get('resid', 'chain %s'%ch))
            len_uniq_res = len(uniq_res)
            if len_uniq_res > len_maj_ch:
                maj_ch = ch
                len_maj_ch = len_uniq_res
        mol_cp.filter('chain %s' % maj_ch)
        return mol_cp
    
    def convert(self):
        os.makedirs("../input/receptor_centered", exist_ok=True)
        if self.conf["docking_type"] == "specified":
            rec_names = np.unique(self.loader.receptor_name)
            file_paths = [os.path.join(self.path, "%s.pdb"%rec_name) for rec_name in rec_names]
        else:
            file_paths =  glob.glob(os.path.join(self.path, "*.pdb"))
        for path in tqdm(file_paths):
            rec_name =  path.split("/")[-1].split(".")[0]
            pdbqt_rec_name =  "../input/receptor_pdbqt/%s.pdbqt"%rec_name
            if not os.path.exists(pdbqt_rec_name):
                try:
                    print("%s is not found, created the new file." % pdbqt_rec_name)
                    self._convert(rec_name)
                except:
                    print("cannot convert into pdbqt: rec_name = %s" % rec_name)
                    pass

class SDFConverter:
    #ここに書く
    def  __init__(self, loader):
        self.loader = loader
    
    def _convert(self, uid):
        obabel = '/home/apps/openbabel/3.1.1/bin/obabel'
        sdf_lig_name = "../input/ligand_dir/%s.sdf" % uid
        pdbqt_lig_name = "../input/ligand_pdbqt/%s.pdbqt" % uid
        cmd = "%s -h -isdf %s -opdbqt -O %s" % (obabel, sdf_lig_name, pdbqt_lig_name)
        os.system(cmd)
        
    def convert(self):
        for uid in tqdm(self.loader.uids): 
            # convert SDF to PDBQT
            pdbqt_lig_name = "../input/ligand_pdbqt/%s.pdbqt" % uid
            if not os.path.exists(pdbqt_lig_name):
                try:
                    print("%s is not found, created the new file." % pdbqt_lig_name)
                    self._convert(uid)
                except:
                    print("cannot convert into pdbqt: sdf_name = %s" %(uid))
                    pass