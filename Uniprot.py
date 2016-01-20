# Projeto_LB
from Bio import SwissProt
from Bio import ExPASy
import os

ids=['O83262', 'O83263', 'O83264', 'O83265', 'O83266', 'O83267', 'O83268', 'O83269', 'O83270', 'O83271', 'O83272', 'O83273', 'ND', 'O83275', 'O83276', 'ND', 'O83278', 'O83279', 'ND', 'O83281', 'O66075', 'O66076', 'O30405', 'O83282', 'O83283', 'O83284', 'O83285', 'ND', 'O83287', 'O83288', 'O83289', 'O83291', 'O83292', 'O83293', 'ND', 'O83295', 'O83296', 'O83297', 'O83298', 'O83299', 'ND', 'O83301', 'ND', 'O83306', 'O83307', 'O83308', 'O83309', 'ND', 'ND', 'ND', 'O83315', 'ND', 'O83317', 'O83318', 'O83319', 'O83320', 'ND', 'O83321', 'O83322', 'O83323', 'O83324', 'ND', 'O83326', 'O83327', 'O83328', 'ND', 'O83330', 'O83331', 'O83332', 'O83334', 'O83335', 'ND', 'ND', 'O83337', 'ND', 'O83341', 'O83342', 'O83343', 'O83344', 'O83345', 'ND', 'O83347', 'O83348', 'ND', 'O83350', 'O83351', 'O83353', 'O83354', 'O83355', 'ND', 'O83357', 'O83358', 'O83359', 'O83360', 'O83361', 'ND', 'O83363', 'O83364', 'O83365', 'O83366', 'O83367', 'ND', 'O83369', 'ND', 'ND', 'O83371', 'O83372', 'O83373', 'ND', 'ND', 'O83377']

class uniprot:
    
    def __init__(self,ids):
        self.ids=ids
        
    def acession(self):
        self.rec=[]
        for ide in self.ids:
            if ide!='ND':
                results=ExPASy.get_sprot_raw(ide)
                rec=SwissProt.read(results)
                self.rec.append(rec)
            else:
                self.rec.append('ND')
        return self.rec
    
    def uniID(self):
        self.uniID=[]
        for rec in self.rec:
            if rec!='ND':
                self.uniID.append(rec.entry_name)
            else:
                self.uniID.append('ND')
        return self.uniID
        
    def Gene_name(self):
        self.gene=[]
        for rec in self.rec:
            if rec!='ND':
                self.gene.append(rec.gene_name)
            else:
                self.gene.append('ND')
        return self.gene
        
    def Revision(self):
        self.rev=[]
        for rec in self.rec:
            if rec!='ND':
                self.rev.append(rec.data_class)
            else:
                self.rev.append('ND')
        return self.rev
    
    def Protein_len(self):
        self.prot_len=[]
        for rec in self.rec:
            if rec!='ND':
                self.prot_len.append(rec.sequence_length)
            else:
                self.prot_len.append('ND')
        return self.prot_len
        
    def Subcel_loc(self):
        self.subcel=[]
        for rec in self.rec:
            if rec!='ND':
                encontrado=False
                com=rec.comments
                for i in range(len(com)):
                    if 'SUBCELLULAR LOCATION:' in com[i]:
                        self.subcel.append(com[i])
                        encontrado=True
                if not encontrado:
                    self.subcel.append('ND')
            else:
                self.subcel.append('ND')
        return self.subcel
        
    def Function(self):
        self.func=[]
        for rec in self.rec:
            if rec!='ND':
                encontrado=False
                com=rec.comments
                for i in range(len(com)):
                    if 'FUNCTION:' in com[i]:
                        self.func.append(com[i])
                        encontrado=True
                if not encontrado:
                    self.func.append('ND')
            else:
                self.func.append('ND')
        return self.func
        
    def Kegg_ort(self):
        self.ko=[]
        for rec in self.rec:
            if rec!='ND':
                enc=False
                ref=rec.cross_references
                for i in range(len(ref)):
                    if ref[i][0]=='KO':
                        self.ko.append(ref[i][1])
                        enc=True
                if not enc:
                    self.ko.append('ND')
            else:
                self.ko.append('ND')
        return self.ko
        
    def Gene_ont(self):
        self.go=[]
        for rec in self.rec:
            if rec!='ND':            
                string=''
                ref=rec.cross_references
                for i in range(len(ref)):
                    if ref[i][0]=='GO':
                        string+=ref[i][1]+'; '
                if string!='':
                    self.go.append(string)
                else:
                    self.go.append('ND')
            else:
                self.go.append('ND')
        return self.go
  
    def GO_names(self):
        self.go_name=[]
        for rec in self.rec:
            if rec!='ND':            
                string=''
                ref=rec.cross_references
                for i in range(len(ref)):
                    if ref[i][0]=='GO':
                        string+=ref[i][2]+'; '
                if string!='':
                    self.go_name.append(string)
                else:
                    self.go_name.append('ND')
            else:
                self.go_name.append('ND')
        return self.go_name
        
    def Patric(self):
        self.pat=[]
        for rec in self.rec:
            if rec!='ND':
                enc=False
                ref=rec.cross_references
                for i in range(len(ref)):
                    if ref[i][0]=='PATRIC':
                        self.pat.append(ref[i][1])
                        enc=True
                if not enc:
                    self.pat.append('ND')
            else:
                self.pat.append('ND')
        return self.pat
        
    def do_protein_table(self):
        uniID=self.uniID()
        gene=self.Gene_name()
        review=self.Revision()
        len_protein=self.Protein_len()
        subcel=self.Subcel_loc()
        func=self.Function()
        kg_ont=self.Kegg_ort()
        g_ont=self.Gene_ont()
        g_ont_names=self.GO_names()
        patric=self.Patric()
        path='protein_table'
        if not os.path.exists(path):
            os.makedirs(path) 
        t="uniID"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(uniID[i])+'\n')     
            temp_file.flush()
            temp_file.close()
        t="gene_name"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(gene[i])+'\n')     
            temp_file.flush()
            temp_file.close()
        t="review"         
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(review[i])+'\n')     
            temp_file.flush()
            temp_file.close() 
        t="len_protein"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(len_protein[i])+'\n')
            temp_file.flush()
            temp_file.close()
        t="subcel"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(subcel[i])+'\n')
            temp_file.flush()
            temp_file.close()
        t="func"
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(func[i])+'\n')
            temp_file.flush()
            temp_file.close()
        t="kg_ont"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(kg_ont[i])+'\n')
            temp_file.flush()
            temp_file.close()
        t="g_ont_names"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(g_ont_names[i])+'\n')
            temp_file.flush()
            temp_file.close()
        t="g_ont"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(g_ont[i])+'\n')
            temp_file.flush()
            temp_file.close()
        t="patric"     
        nome=t+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(uniID)):
                temp_file.writelines(str(patric[i])+'\n')
            temp_file.flush()
            temp_file.close()
