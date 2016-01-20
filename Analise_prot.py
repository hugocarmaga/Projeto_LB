from Bio import SeqIO
from Bio import Entrez
import os


class Analise_prot:
    
    def __init__(self, Lista_GI, database='protein'):
        self.acession=Lista_GI
        self.db=database
        self.records=None
    
    def get_acession(self):
        self.records=[]
        for i in self.acession:
            if i != 'ND':
                handle = Entrez.efetch(db=self.db, rettype="gp", retmode='text', id=i)
                record=SeqIO.read(handle, 'genbank')
                self.records.append(record)
                handle.close()
            else:
                self.records.append('ND')
        return self.records
    
    def obtain_CDD(self):
        self.cdds=[]
        self.get_acession()
        for record in self.records:
            if record != 'ND':
                CDDS=''
                for feat in record.features:
                    if feat.type=='Region':
                        if 'db_xref' in feat.qualifiers:
                            CDDS=CDDS + '   ' + str(feat.qualifiers['db_xref'])
                self.cdds.append(CDDS)
            else:
                self.cdds.append('ND')
        return self.cdds
    
    def obtain_notes(self):
        self.note_cdds=[]
        self.get_acession()
        for record in self.records:
            if record != 'ND':
                notee=''
                for feat in record.features:
                    if feat.type=='Region':
                        if 'note' in feat.qualifiers:
                            notee=notee + '|' + '   ' + str(feat.qualifiers['note'])
                self.note_cdds.append(notee)
            else:
                self.note_cdds.append('ND')
        return self.note_cdds
        
    
    def do_CCD_table(self):
        cdds=self.obtain_CDD()
        path='cdd_table'   
        if not os.path.exists(path):
            os.makedirs(path)
        nome=path+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(cdds)):
                temp_file.writelines(str(cdds[i])+'\n')
            temp_file.flush()
            temp_file.close()
    
    def do_notes_table(self):
        notes_cdds=self.obtain_notes()
        path='notes_table'   
        if not os.path.exists(path):
            os.makedirs(path)
        nome=path+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(notes_cdds)):
                temp_file.writelines(str(notes_cdds[i])+'\n')
            temp_file.flush()
            temp_file.close()
