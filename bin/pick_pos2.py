# -*- coding=utf-8 -*-
# 2021.0706
# @zlk 
# zhanglk960127@163.com
# gffread 提取cds

import os


def get_pep(gffread,out_file):
    os.system("%s %s/species.gff -g %s/species.fa -x %s/species.gff.cds" %(gffread,out_file,out_file,out_file))
    in_cds_file=open(out_file+'/species.gff.cds','r')
    out_pep_file=open(out_file+'/species.gff.pep','w')
    cds_dict={}
    for line in in_cds_file:
        if line.startswith('>'):
            cds_gene_id=line.strip()
            cds_dict[cds_gene_id]=''
        else:
            cds_dict[cds_gene_id]+=line.strip()
    for key,value in cds_dict.items():
        out_pep_file.write(key+'\n'+rna_translate(value)+'\n')
    

        


