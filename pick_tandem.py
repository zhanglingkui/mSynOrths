# -*- coding=utf-8 -*-
# 2020.1008
# @zlk
# 去除tandem基因，保留cds最长的一个基因
import os
import multiprocessing
def filter_blast(blast_file):
    fiter_file=open( blast_file,'r')
    fiter_gene_dict={}
    fiter_gene_list=[]
    for fiter_line in fiter_file:
        fiter_line_list=fiter_line.strip().split('\t')
        if fiter_line_list[0] in fiter_gene_dict.keys():
            if fiter_line_list[1] in fiter_gene_dict[fiter_line_list[0]]:
                continue
            else:
                fiter_gene_list.append(fiter_line)
            fiter_gene_dict[fiter_line_list[0]].append(fiter_line_list[1])
        else:
            fiter_gene_dict[fiter_line_list[0]]=[fiter_line_list[1]]
            fiter_gene_list.append(fiter_line)
    fiter_file.close()
    fiter_file2=open(blast_file,'w')
    for fiter_line in fiter_gene_list:
        fiter_file2.write(fiter_line)
    fiter_file2.close()
    del fiter_gene_dict
    del fiter_gene_list

def do_blast(query_prot,ref_prot,work_path,tools,threads_num,evalue,out_name): 
    if tools == 'diamond':
                commond11 = 'mkdir ' + work_path + '/db/'
                os.system(commond11)
                commond1 = 'diamond makedb --in ' + ref_prot + ' -d ' + work_path + '/db/ref'
                print(commond1)
                os.system(commond1)
                commond2 = 'diamond blastp -d ' + work_path + '/db/ref -q ' + query_prot + ' --quiet --sensitive -e '+evalue+'  -p ' + str(
                    threads_num) + ' -o ' + work_path + out_name
                os.system(commond2)
    else:
        commond1 = 'makeblastdb -in ' + ref_prot + ' -dbtype prot -parse_seqids -out ' + work_path + '/db/ref'
        os.system(commond1)
        commond2 = 'blastp -db ' + work_path + '/db/ref -query ' + query_prot + '  -num_threads ' + str(
            threads_num) + ' -outfmt 6 -evalue '+evalue+' -out ' + work_path + out_name
        os.system(commond2)
    filter_blast(work_path + out_name)  

def self_blast(out_folder,tools,threads_num):
    for mSynF in os.listdir(out_folder):
        if mSynF.startswith('mSynF'):
            work_path=out_folder+'/'+mSynF+'/'
            if 'self_blast' in os.listdir(work_path):
                continue
            prot_file=work_path+'/species.gff_pep'
            do_blast(prot_file,prot_file,work_path,tools,threads_num,'1e-5','self_blast')

def pick_tandem(pos_file,self_blast_file,tandem_file,prot_file):
    tandem_num=5
    op_pos_file=open(pos_file,'r')
    gene_pos_length_dict={}
    pos_index=0
    for line in op_pos_file:
        pos_index+=1
        line_list=line.strip().split('\t')
        gene_pos_length_dict[line_list[0]]=[pos_index,int(line_list[2])/3]
    op_pos_file.close()
    op_blast=open(self_blast_file,'r')
    tandem_dict={}
    for blast_line in op_blast:
        blast_line_list=blast_line.strip().split('\t')
        try:
            if float(blast_line_list[2])<30:
                continue
            elif int(blast_line_list[3])/gene_pos_length_dict[blast_line_list[0]][1]<0.3 or int(blast_line_list[3])/gene_pos_length_dict[blast_line_list[1]][1]<0.3:
                continue
            elif 0<gene_pos_length_dict[blast_line_list[1]][0] - gene_pos_length_dict[blast_line_list[0]][0] <=tandem_num:
                if blast_line_list[0] in tandem_dict.keys():
                    tandem_dict[blast_line_list[0]].append(blast_line_list[1])
                else:
                    tandem_dict[blast_line_list[0]]=[blast_line_list[1]]
        except KeyError:
            continue
    op_blast.close()
    # del gene_pos_length_dict
    del_list=[]
    tandem_dict_new={}
    for key,value in tandem_dict.items():
        if key in del_list:
            continue
        for i in value:
            if (i not in del_list) and (i in tandem_dict.keys()):
                for j in tandem_dict[i]:
                    if j not in value:
                        value.append(j)
                # value.append(tandem_dict[i])
                del_list.append(i)
                # value=list(set(value))
        tandem_dict_new[key]=value
    del tandem_dict
    
    op_tandem_file=open(tandem_file,'w')
    tandem_list=[]
    for key,value in tandem_dict_new.items():
        aplist=[]
        aplist.append(key)
        for i in value:
            aplist.append(i)
        tandem_list.append(aplist)
    del_tandem_list=[]
    for i_list in tandem_list:
        for i in sorted(i_list,key=lambda x: gene_pos_length_dict[x][1])[1:]:
            del_tandem_list.append(i)
        for i in sorted(i_list,key=lambda x: gene_pos_length_dict[x][1]):
            op_tandem_file.write(i+'\t')
        op_tandem_file.write('\n')
    op_tandem_file.close()
    op_pos_file=open(pos_file,'r')
    wr_pos_file=open(pos_file+'_DelTandem','w')
    for line in op_pos_file:
        line_list=line.strip().split('\t')
        if line_list[0] in del_tandem_list:
            continue
        else:
            wr_pos_file.write(line)
    op_pos_file.close()
    wr_pos_file.close()
    op_pep_file=open(prot_file,'r')
    wr_pep_file=open(prot_file+'_DelTandem','w')
    flag=1
    for line in op_pep_file:
        if line.startswith('>'):
            if line[1:].strip() in del_tandem_list:
                flag=0
            elif line[1:].strip() not in gene_pos_length_dict:
                flag=0
            else:
                flag=1
        if flag==1:
            wr_pep_file.write(line)
    op_pep_file.close()
    wr_pep_file.close()
    del tandem_dict_new
    del gene_pos_length_dict
    del tandem_list
def run_pick_tandem(out_folder,threads_num):
    pool = multiprocessing.Pool(threads_num)
    for mSynF in os.listdir(out_folder):
        if mSynF.startswith('mSynF'):
            work_path=out_folder+'/'+mSynF+'/'
            print(work_path)
            if 'tandem_array.txt' in os.listdir(work_path):
                result=pool.apply_async(pick_tandem,(work_path+'species.gff_pos',work_path+'/self_blast',work_path+'/tandem_array.txt',work_path+'/species.gff_pep'))
                # continue
            else:
                result=pool.apply_async(pick_tandem,(work_path+'species.gff_pos',work_path+'/self_blast',work_path+'/tandem_array.txt',work_path+'/species.gff_pep'))
    result.get()
    pool.close()
    pool.join()    

if __name__ == "__main__":
    import sys
    # self_blast(sys.argv[1],sys.argv[2],10)
    # pick_tandem(sys.argv[1],sys.argv[2],sys.argv[3])
    run_pick_tandem(sys.argv[1],int(sys.argv[2]))