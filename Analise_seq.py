from Bio import SeqIO
from Bio import Entrez
import re
import os
import urllib

class Analise_Seq:


    def __init__(self, NC_number, database='nucleotide', start_position=None, end_position=None):
        self.acession=NC_number
        self.db=database
        self.start=start_position
        self.end=end_position
        self.record=None
    
    def get_acession(self):
        handle = Entrez.efetch(self.db, rettype='gb', retmode='text', seq_start=self.start, seq_stop=self.end, id=self.acession)
        self.record=SeqIO.read(handle, 'genbank')
        handle.close()
        return self.record
    
    def write_record(self, filename):
        if self.record != None:
            SeqIO.write(self.record, filename, 'genbank')
        else:
            self.get_acession()
            SeqIO.write(self.record, filename, 'genbank')
    
    def print_record(self):
        if self.record != None:
            print(self.record) 
        else:
            self.get_acession()
            print(self.record)
            
    def print_features(self):
        if self.record != None:
            for feat in self.record.features:
                print (feat.qualifiers) 
        else:
            self.get_acession()
            for feat in self.record.features:
                print (feat.qualifiers)
                
    def escreve_fasta(self):
        sequencias=self._get_seq()
        locus_tag=self._get_locustag()
        gene_ids=self._get_geneid()
        for i in range(len(gene_ids)): 
            path='FASTA_BLAST'   
            if not os.path.exists(path):
                os.makedirs(path)
            nome=str(gene_ids[i])+".fasta"
            with open(os.path.join(path, nome), 'w') as temp_file:
                temp_file.writelines('>' +gene_ids[i]+ '  ' + locus_tag[i]+'\n')
                temp_file.writelines(sequencias[i])
                temp_file.flush()
                temp_file.close()
    
    def do_protein_table(self):
        locus_tag=self._get_locustag()
        gene_ids=self._get_geneid()
        protein_ids=self._get_proteinid()
        gis=self._get_gi()
        locations=self._get_local()
        ids=self.uniprot_ID()
        ec=self.ec_number()
        #ids=['O83262', 'O83263', 'O83264', 'O83265', 'O83266', 'O83267', 'O83268', 'O83269', 'O83270', 'O83271', 'O83272', 'O83273', 'ND', 'O83275', 'O83276', 'ND', 'O83278', 'O83279', 'ND', 'O83281', 'O66075', 'O66076', 'O30405', 'O83282', 'O83283', 'O83284', 'O83285', 'ND', 'O83287', 'O83288', 'O83289', 'O83291', 'O83292', 'O83293', 'ND', 'O83295', 'O83296', 'O83297', 'O83298', 'O83299', 'ND', 'O83301', 'ND', 'O83306', 'O83307', 'O83308', 'O83309', 'ND', 'ND', 'ND', 'O83315', 'ND', 'O83317', 'O83318', 'O83319', 'O83320', 'ND', 'O83321', 'O83322', 'O83323', 'O83324', 'ND', 'O83326', 'O83327', 'O83328', 'ND', 'O83330', 'O83331', 'O83332', 'O83334', 'O83335', 'ND', 'ND', 'O83337', 'ND', 'O83341', 'O83342', 'O83343', 'O83344', 'O83345', 'ND', 'O83347', 'O83348', 'ND', 'O83350', 'O83351', 'O83353', 'O83354', 'O83355', 'ND', 'O83357', 'O83358', 'O83359', 'O83360', 'O83361', 'ND', 'O83363', 'O83364', 'O83365', 'O83366', 'O83367', 'ND', 'O83369', 'ND', 'ND', 'O83371', 'O83372', 'O83373', 'ND', 'ND', 'O83377']
        path='protein_table'   
        if not os.path.exists(path):
            os.makedirs(path)
        nome=path+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(gene_ids)):
                temp_file.writelines('chr    NC_000919.1'+ '    ' +locations[i]+ '    '+gene_ids[i]+ '    '+locus_tag[i]+ '    '+protein_ids[i]+ '    '+gis[i]+ '     '+ids[i]+'     '+ec[i]+'\n')
            temp_file.flush()
            temp_file.close()
    
    def ec_number(self):
        self.ecNumbs=[]
        for i in range(len(self.record.features)):
            my_gene = self.record.features[i]
            if my_gene.type == "CDS": 
                if "EC_number" in my_gene.qualifiers:
                    self.ecNumbs.append(my_gene.qualifiers["EC_number"][0])
                else:
                    self.ecNumbs.append('ND')
        return self.ecNumbs
    
    def uniprot_ID(self):
        self.IDs=[]
        protein_ids=self._get_proteinid()
        for j in range(len(protein_ids)):#get uniprotAC for each proteinID from url
        #the html page has in a consistent position the uniprotAC
            if protein_ids[j] != 'ND':
                data = urllib.request.urlopen("http://www.uniprot.org/uniprot/?query="+protein_ids[j]+"&sort=score").read()
                m=re.search('id="O', str(data))
                if m != None:
                    ini=m.span()[1]-1
                    fim=m.span()[1]+5
                    self.IDs.append(str(data)[ini:fim])
                else:
                    self.IDs.append('ND')
            else:
                self.IDs.append('ND')
        #ids=['O83262', 'O83263', 'O83264', 'O83265', 'O83266', 'O83267', 'O83268', 'O83269', 'O83270', 'O83271', 'O83272', 'O83273', 'ND', 'O83275', 'O83276', 'ND', 'O83278', 'O83279', 'ND', 'O83281', 'O66075', 'O66076', 'O30405', 'O83282', 'O83283', 'O83284', 'O83285', 'ND', 'O83287', 'O83288', 'O83289', 'O83291', 'O83292', 'O83293', 'ND', 'O83295', 'O83296', 'O83297', 'O83298', 'O83299', 'ND', 'O83301', 'ND', 'O83306', 'O83307', 'O83308', 'O83309', 'ND', 'ND', 'ND', 'O83315', 'ND', 'O83317', 'O83318', 'O83319', 'O83320', 'ND', 'O83321', 'O83322', 'O83323', 'O83324', 'ND', 'O83326', 'O83327', 'O83328', 'ND', 'O83330', 'O83331', 'O83332', 'O83334', 'O83335', 'ND', 'ND', 'O83337', 'ND', 'O83341', 'O83342', 'O83343', 'O83344', 'O83345', 'ND', 'O83347', 'O83348', 'ND', 'O83350', 'O83351', 'O83353', 'O83354', 'O83355', 'ND', 'O83357', 'O83358', 'O83359', 'O83360', 'O83361', 'ND', 'O83363', 'O83364', 'O83365', 'O83366', 'O83367', 'ND', 'O83369', 'ND', 'ND', 'O83371', 'O83372', 'O83373', 'ND', 'ND', 'O83377']
        return self.IDs
    
    def _get_seq(self, tipo='CDS'):
        seq=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:
                seq.append(self.record.features[i].extract(self.record.seq))
        return seq
    
    def _get_local(self,tipo='CDS'):
        local_list=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:
                local=str(self.record.features[i].location)
                local_right=''
                for char in local:
                    if char in '1234567890+-':
                        local_right=local_right+char
                    if char in ':(':
                        local_right=local_right+'   '
                local_list.append(local_right)
        return local_list
    
    def _get_geneid(self,tipo='CDS'):
        geneid=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:        
                gene_id=''
                if re.search('GI.[0-9]*',self.record.features[i].qualifiers["db_xref"][0]) != None:
                    gene_id=self.record.features[i].qualifiers["db_xref"][1]
                else:
                    gene_id=self.record.features[i].qualifiers["db_xref"][0]
                gene_id_rigth=''
                for cha in gene_id:
                    if cha in '1234567890':
                        gene_id_rigth=gene_id_rigth+cha
                geneid.append(gene_id_rigth)
        return geneid
    def _get_gi(self,tipo='CDS'):
        gi=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:        
                gene_gi=''
                if re.search('GI.[0-9]*',self.record.features[i].qualifiers["db_xref"][0]) != None:
                    gene_gi=self.record.features[i].qualifiers["db_xref"][0]
                    gene_gi_rigth=''
                    for cha in gene_gi:
                        if cha in '1234567890':
                            gene_gi_rigth=gene_gi_rigth+cha
                    gi.append(gene_gi_rigth)
                else:
                    gi.append('ND')
        return gi
        
    def _get_locustag(self,tipo='CDS'):
        locus=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:
                locus.append(self.record.features[i].qualifiers["locus_tag"][0])
        return locus
    
    def _get_proteinid(self,tipo='CDS'):
        proteinid=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:
                if self.record.features[i].qualifiers.get("protein_id"):
                    proteinid.append(self.record.features[i].qualifiers["protein_id"][0])
                else:
                    proteinid.append('ND')
        return proteinid
    
    def _get_function(self,tipo='CDS'):
        if self.record == None:
            self.get_acession()
        self.functions=[]
        for i in range(len(self.record.features)):
            my_gene = self.record.features[i]
            if my_gene.type == tipo: 
                if "product" in my_gene.qualifiers:
                    self.functions.append(my_gene.qualifiers["product"][0])
                else:
                    self.functions.append('ND')
        return self.functions
    
    def _get_gis_hipothetical_P(self,tipo='CDS'):
        gi_p=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo:        
                gene_gi=''
                if re.search('GI.[0-9]*',self.record.features[i].qualifiers["db_xref"][0]) != None:
                    if 'hypothetical protein' in self.record.features[i].qualifiers["product"]:
                        gene_gi=self.record.features[i].qualifiers["db_xref"][0]
                        gene_gi_rigth=''
                        for cha in gene_gi:
                            if cha in '1234567890':
                                gene_gi_rigth=gene_gi_rigth+cha
                        gi_p.append(gene_gi_rigth)
        return gi_p
    
    def _get_geneid_hipothetical_P(self,tipo='CDS'):
        geneid_p=[]
        if self.record == None:
            self.get_acession()
        for i in range(len(self.record.features)):
            if self.record.features[i].type == tipo: 
                if 'hypothetical protein' in self.record.features[i].qualifiers["product"]:       
                    gene_id=''
                    if re.search('GI.[0-9]*',self.record.features[i].qualifiers["db_xref"][0]) == None:
                        gene_id=self.record.features[i].qualifiers["db_xref"][0]
                        gene_id_rigth=''
                        for cha in gene_id:
                            if cha in '1234567890':
                                gene_id_rigth=gene_id_rigth+cha
                        geneid_p.append(gene_id_rigth)
        return geneid_p
        
    def do_function_table(self):
        functions=self._get_function()
        path='fucntion_table'   
        if not os.path.exists(path):
            os.makedirs(path)
        nome=path+".txt"
        with open(os.path.join(path, nome), 'w') as temp_file:
            for i in range(len(functions)):
                temp_file.writelines(functions[i]+ '\n')
            temp_file.flush()
            temp_file.close()
