import os
    

class AFDBFetcher:
    def __init__(self, loader, afdb_path):
        self.loader =  loader
        self.afdb_path =  afdb_path

    def _fetch(self, uniprot_id):
        filename = "AF-%s-F1-model_v4.pdb.gz" % uniprot_id 
        source_file = os.path.join(self.afdb_path, filename)
        os.system("gzip -c -d %s > %s"%(source_file, "../input/receptor_dir/%s"%(uniprot_id+".pdb")))
        
    def fetch(self):
        for uniprot_id in self.loader.uniprot_ids:
            self._fetch(uniprot_id)

