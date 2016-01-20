from Analise_seq import Analise_Seq
from Analise_prot import Analise_prot
from Uniprot import uniprot
from blast import blast_p
from blast import blast_n

def main_analise_lit():
    treponema=search_pubmed('Treponema pallidum pallidum genome')
    treponema.records_number()
    treponema.write_records('teste.txt')
    sifilis=search_pubmed('Treponema pallidum pallidum sifilis')
    sifilis.records_number()
    sifilis.write_records('test_sifilis.txt')

def main_analise_seq():
    an_seq=Analise_Seq('NC_000919', start_position=241801, end_position=381900)
    an_seq.write_record('NC_000919.gb')
    an_seq.do_protein_table()
    an_seq.do_function_table()
    gi=an_seq._get_gi()
    print(gi)
    products=an_seq._get_function()
    print(products)
    an_seq.do_protein_table()
    an_seq.escreve_fasta()
    lista_prots=an_seq._get_proteinid()
    print(lista_prots)
    u_ID=an_seq.uniprot_ID()
    print(u_ID)
    ids=['O83262', 'O83263', 'O83264', 'O83265', 'O83266', 'O83267', 'O83268', 'O83269', 'O83270', 'O83271', 'O83272', 'O83273', 'ND', 'O83275', 'O83276', 'ND', 'O83278', 'O83279', 'ND', 'O83281', 'O66075', 'O66076', 'O30405', 'O83282', 'O83283', 'O83284', 'O83285', 'ND', 'O83287', 'O83288', 'O83289', 'O83291', 'O83292', 'O83293', 'ND', 'O83295', 'O83296', 'O83297', 'O83298', 'O83299', 'ND', 'O83301', 'ND', 'O83306', 'O83307', 'O83308', 'O83309', 'ND', 'ND', 'ND', 'O83315', 'ND', 'O83317', 'O83318', 'O83319', 'O83320', 'ND', 'O83321', 'O83322', 'O83323', 'O83324', 'ND', 'O83326', 'O83327', 'O83328', 'ND', 'O83330', 'O83331', 'O83332', 'O83334', 'O83335', 'ND', 'ND', 'O83337', 'ND', 'O83341', 'O83342', 'O83343', 'O83344', 'O83345', 'ND', 'O83347', 'O83348', 'ND', 'O83350', 'O83351', 'O83353', 'O83354', 'O83355', 'ND', 'O83357', 'O83358', 'O83359', 'O83360', 'O83361', 'ND', 'O83363', 'O83364', 'O83365', 'O83366', 'O83367', 'ND', 'O83369', 'ND', 'ND', 'O83371', 'O83372', 'O83373', 'ND', 'ND', 'O83377']
    an_gene=Analise_prot(ids)
    an_gene.get_acession()
    an_gene.write_records()

def main_analise_prot():
    an_seq=Analise_Seq('NC_000919', start_position=241801, end_position=381900)
    gis=an_seq._get_gi()
    an_prot=Analise_prot(gis)
    an_prot.do_CCD_table()
    an_prot.do_notes_table()
    
def main_analise_uniprot():
    an_seq=Analise_Seq('NC_000919', start_position=241801, end_position=381900)
    u_ID=an_seq.uniprot_ID()
    an_uniprot=uniprot(u_ID)
    print(an_uniprot)

def main_blast_p():
    an_seq=Analise_Seq('NC_000919', start_position=241801, end_position=381900)
    gis=an_seq._get_gis_hipothetical_P()
    #retorna gis que pertencam a proteinas hipoteticas
    print(len(gis))
    print(gis)
    blast_p(gis)

def main_blast_n():
    an_seq=Analise_Seq('NC_000919', start_position=241801, end_position=381900)
    gene_ids=an_seq._get_geneid_hipothetical_P()
    #retorna gene_ids que pertencem a proteinas hipoteticas e nao tem gis
    print(len(gene_ids))
    print(gene_ids)
    blast_n(gene_ids)
