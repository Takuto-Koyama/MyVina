import os
import json 

import glob
                 

class PairwiseDocking:
    def __init__(self, config, loader):
        self.config = config
        self.loader = loader
    
    def _dock(self, ligand_name, receptor_name, conf, bconf):
        qvina = '/home/takutoyo/FUJITSU/qvina/bin/qvina-w'
        rec_pdbqt = "../input/receptor_pdbqt/%s.pdbqt"%receptor_name
        lig_pdbqt = "../input/ligand_pdbqt/%s.pdbqt"%ligand_name
        if not os.path.exists(os.path.join(conf["output_dir"],receptor_name,"%s.pdbqt"%ligand_name)):
            try:        
                os.system("%s --receptor %s --ligand %s --seed 0\
                            --center_x %s --center_y %s --center_z %s\
                            --size_x %s --size_y %s --size_z %s\
                            --out %s --num_modes %s --exhaustiveness %s"%\
                            (qvina, rec_pdbqt, lig_pdbqt, 
                            bconf["center_x"], bconf["center_y"], bconf["center_z"], 
                            bconf["size_x"], bconf["size_y"], bconf["size_z"],
                            os.path.join(conf["output_dir"],receptor_name,"%s.pdbqt"%(ligand_name)), 
                            conf["num_modes"], conf["exhaustiveness"]))
            except:
                print("Counld not execute Quick Vina-W for receptor_name = %s and ligand_name = %s "%(receptor_name, ligand_name))
        else:
            print("%s-%s.pdbqt already exists"%(receptor_name, ligand_name))
    
    def dock(self):
        if self.config["ligand_type"] == "smiles":
            if self.config["protein_type"] == "sequence":
                # csv loader
                raise NotImplementedError
            elif self.config["protein_type"] == "pdb":
                # liglist loader
                raise NotImplementedError
            else:
                raise ValueError("protein_type must be sequence, pdb, or pdbqt")
        elif self.config["ligand_type"] == "sdf":
            # liglist loader
            receptor_paths = glob.glob("../input/receptor_pdbqt/*.pdbqt")
            
            for receptor_path in receptor_paths:
                receptor_name = receptor_path.split("/")[-1].split(".")[0]
                os.makedirs("../output/%s"%receptor_name, exist_ok=True)
                bconf = "../input/receptor_pdbqt/%s_box.json"%receptor_name
                with open(bconf, "r") as f:
                    bconf = json.load(f)
                for ligand_name in self.loader.uids:
                    self._dock(ligand_name, receptor_name, self.config, bconf)
        
        
class SpecifiedDocking:
    def __init__(self, config, loader):
        self.config = config
        self.loader = loader
        
    def _dock(self, ligand_name, receptor_name, conf, bconf):
        qvina = '/home/takutoyo/FUJITSU/qvina/bin/qvina-w'
        rec_pdbqt = "../input/receptor_pdbqt/%s.pdbqt"%receptor_name
        lig_pdbqt = "../input/ligand_pdbqt/%s.pdbqt"%ligand_name
        if not os.path.exists(os.path.join(conf["output_dir"],receptor_name,"%s.pdbqt"%ligand_name)):
            try:
                print(os.path.join(conf["output_dir"],receptor_name,"%s.pdbqt"%(ligand_name)))        
                os.system("%s --receptor %s --ligand %s --seed 0\
                            --center_x %s --center_y %s --center_z %s\
                            --size_x %s --size_y %s --size_z %s\
                            --out %s --num_modes %s --exhaustiveness %s"%\
                            (qvina, rec_pdbqt, lig_pdbqt, 
                            bconf["center_x"], bconf["center_y"], bconf["center_z"], 
                            bconf["size_x"], bconf["size_y"], bconf["size_z"],
                            os.path.join(conf["output_dir"],receptor_name,"%s.pdbqt"%(ligand_name)), 
                            conf["num_modes"], conf["exhaustiveness"]))
            except:
                print("Counld not execute Quick Vina-W for receptor_name = %s and ligand_name = %s "%(receptor_name, ligand_name))
        else:
            print("%s-%s.pdbqt already exists"%(receptor_name, ligand_name))
    
    def dock(self):
        if self.config["ligand_type"] == "smiles":
            if self.config["protein_type"] == "sequence":
                for uniprot_id in self.loader.uniprot_ids:
                    df_uniprot = self.loader.df[self.loader.df["uniprot_id"]==uniprot_id]
                    for drug_id in df_uniprot["drug_id"]:
                        receptor_name = uniprot_id
                        ligand_name = drug_id
                        self._dock(ligand_name, receptor_name, self.config)
            else:
                raise NotImplementedError
        
        elif self.config["ligand_type"] == "sdf":
            if self.config["protein_type"] == "pdb":
                for ligand_name, receptor_name in zip(self.loader.ligand_name, self.loader.receptor_name):
                    os.makedirs("../output/%s"%receptor_name, exist_ok=True)
                    bconf = "../input/receptor_pdbqt/%s_box.json"%receptor_name
                    try:
                        with open(bconf, "r") as f:
                            bconf = json.load(f)
                        self._dock(ligand_name, receptor_name, self.config, bconf)
                    except:
                        pass
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError