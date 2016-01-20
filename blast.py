from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def blast_p(gis):
    for ide in gis:
        if ide != 'ND':
            nome=str(ide)+".xml"
            save_file=open(nome, 'w')
            result=NCBIWWW.qblast("blastp","swissprot",str(ide),hitlist_size=5)
            save_file.write(result.read())
            result.close()
            save_file.close()
            
            result_file=open(nome, 'r')
            blast= NCBIXML.parse(result_file)

                
            E_VALUE_THRESH = 0.04
            for blast_record in blast:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            print ('\n')
                            print ('****Alignment****')
                            print("\n")
                            print("Acession numbers: " ,alignment.accession)
                            print("ID do alinhamento: " ,alignment.hit_id)
                            print("Descrição do alinhamento: " ,alignment.hit_def)
                            print ('sequence:', alignment.title)
                            print('length:', alignment.length)
                            print('e value:', hsp.expect)
                            print ('e value:', hsp.expect)
                            print (hsp.query[0:75])
                            print (hsp.match[0:75])
                            print (hsp.sbjct[0:75])
        else:
            print ('\n')
            print ('****Alignment****')
            print("\n")
            print('ND')

def blast_n(gene_id):
    for ide in gene_id:
        if ide != 'ND':
            nome=str(ide)+".xml"
            save_file=open(nome, 'w')
            result=NCBIWWW.qblast("blastn","nt",str(ide),hitlist_size=5)
            save_file.write(result.read())
            result.close()
            save_file.close()
            
            result_file=open(nome, 'r')
            blast= NCBIXML.parse(result_file)

                
            E_VALUE_THRESH = 0.04
            for blast_record in blast:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            print ('\n')
                            print ('****Alignment****')
                            print("\n")
                            print("Acession numbers: " ,alignment.accession)
                            print("ID do alinhamento: " ,alignment.hit_id)
                            print("Descrição do alinhamento: " ,alignment.hit_def)
                            print ('sequence:', alignment.title)
                            print('length:', alignment.length)
                            print('e value:', hsp.expect)
                            print ('e value:', hsp.expect)
                            print (hsp.query[0:75])
                            print (hsp.match[0:75])
                            print (hsp.sbjct[0:75])
        else:
            print ('\n')
            print ('****Alignment****')
            print("\n")
            print('ND')
