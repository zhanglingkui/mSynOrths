# -*- coding=utf-8 -*-
# 2020.1007
# @zlk
# 提取蛋白质序列，和posfile
import sys
import re
import os
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
# 在gff文件第九列提取出其ID


def pick_mRNA(gff_file,fa_file,out_folder):
    
    def get_pep(gffread,out_file):
        os.system("%s %s/species.gff -g %s/species.fa -S -y %s/species.gff_pep" %(gffread,out_file,out_file,out_file))
    get_pep('gffread',out_folder)
    def pick_id(line9):
        list9 = line9.split(';')
        for i in list9:
            if i.startswith('ID='):
                if ',' in i:
                    return i.split(',')[0][3:]
                else:
                    return i[3:]

    def pick_parent(line9):
        list9 = line9.split(';')
        for i in list9:
            if i.startswith('Parent='):
                if ',' in i:
                    return i.split(',')[0][7:]
                else:
                    return i[7:]
    op_gff_file=open(gff_file,'r')
    gff_list=[]
    for gff_line in op_gff_file:
        gff_line_list=gff_line.strip().split('\t')
        if len(gff_line_list)==9:
            if gff_line_list[2]=='gene' or gff_line_list[2]=='CDS' or gff_line_list[2]=='mRNA':
                gff_list.append(gff_line_list)
    op_gff_file.close()
    gene_index=0
    mRNA_index=0
    for gff_i_list in gff_list[:200]:
        if gff_i_list[2]=='gene':
            gene_index=1
        elif gff_i_list[2]== 'mRNA':
            mRNA_index=1
    if gene_index:
        if mRNA_index:
            gff_type='gene_mRNA'
        else:
            gff_type='gene'
    else:
        if mRNA_index:
            gff_type='mRNA'
        else:
            print(gff_file+ ' format is wrong')
            sys.exit()
    if gff_type=='mRNA' or gff_type=='gene':
        gene_dict={}
        for gff_line in gff_list:
            if gff_line[2]=='CDS':
                # gff_id=pick_id(gff_line[8])
                gff_parent=pick_parent(gff_line[8])
                if int(gff_line[4])> int(gff_line[3]):
                    aplist=[int(gff_line[3]),int(gff_line[4]),gff_line[0]]
                else:
                    aplist=[int(gff_line[4]),int(gff_line[3]),gff_line[0]]
                if gff_parent in gene_dict.keys():
                    gene_dict[gff_parent][1]+=aplist[1]-aplist[0]+1
                    gene_dict[gff_parent].append(aplist)
                else:
                    gene_dict[gff_parent]=[gff_line[6],aplist[1]-aplist[0]+1]
                    gene_dict[gff_parent].append(aplist)
    else:
        gene_dict={}
        gene_mrna_dict={}
        for gff_line in gff_list:
            if gff_line[2]=='CDS':
                # gff_id=pick_id(gff_line[8])
                gff_parent=pick_parent(gff_line[8])
                if int(gff_line[4])> int(gff_line[3]):
                    aplist=[int(gff_line[3]),int(gff_line[4]),gff_line[0]]
                else:
                    aplist=[int(gff_line[4]),int(gff_line[3]),gff_line[0]]
                if gff_parent in gene_dict.keys():
                    gene_dict[gff_parent][1]+=aplist[1]-aplist[0]+1
                    gene_dict[gff_parent].append(aplist)
                else:
                    gene_dict[gff_parent]=[gff_line[6],aplist[1]-aplist[0]+1]
                    gene_dict[gff_parent].append(aplist)
            elif gff_line[2]=='mRNA':
                gff_parent=pick_parent(gff_line[8])
                ##############
                gff_id=pick_id(gff_line[8])
                if gff_parent:
                    a=1
                else:
                   gff_parent= gff_id
                if gff_parent in gene_mrna_dict.keys():
                    gene_mrna_dict[gff_parent].append(gff_id)
                else:
                    gene_mrna_dict[gff_parent]=[gff_id]
        gene_dict_new={}
        for key,value in gene_mrna_dict.items():
            if len(value)==1:
                gene_dict_new[value[0]]=gene_dict[value[0]] 
            else:
                value_sort=sorted(value ,key=lambda x:gene_dict[x][1])
                gene_dict_new[value_sort[-1]]=gene_dict[value_sort[-1]] 
        gene_dict=gene_dict_new
    gff_end=gff_file.split('/')[-1]
    outfile_pos=open(out_folder+'/'+gff_end+'_pos','w')
    chr_gene_dict={}
    for key,value in gene_dict.items():
        # fa_seq=''
        # for cds in sorted(value[2:],key=lambda x:x[0]):
        #     x=0
        cds=value[-1]
        if cds[2] in chr_gene_dict.keys():
            chr_gene_dict[cds[2]].append([key]+value)
        else:
            chr_gene_dict[cds[2]]=[[key]+value]
    for key,value in chr_gene_dict.items():
        ## sort gff###
        for i in sorted(value,key=lambda x:x[3][0]):
            for j in i[:3]:
                outfile_pos.write(str(j)+'\t')
            outfile_pos.write(key+'\t')
            for j in i[3:]:
                outfile_pos.write(str(j[0])+':'+str(j[1])+'\t')
            outfile_pos.write('\n')


def run_pick_pep(threads,out_folder):
    pool = multiprocessing.Pool(threads)
     
    for i in os.listdir(out_folder):
        if i.startswith('mSynF'):
            i_path=out_folder+'/'+i
            # if 'species.gff_pep' in os.listdir(i_path):
            #     continue
            result = pool.apply_async(pick_mRNA,(i_path+'/species.gff',i_path+'/species.fa',i_path))
    result.get()
    pool.close()
    pool.join()


if __name__ == "__main__":
    # pick_mRNA(sys.argv[1],sys.argv[2],sys.argv[3])
    run_pick_pep(int(sys.argv[1]),sys.argv[2])



    
    

    


