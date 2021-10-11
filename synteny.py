# -*- coding=utf-8 -*-
# 2020.1009
# @zlk
import os
import pick_tandem
import syntenic_segment
import multiprocessing

def together_blast(out_folder,tools,threads_num,evalue):
    species_num=0
    for i in os.listdir(out_folder):
        if i.startswith('mSynF'):
            species_num+=1
    for x in range(1,species_num,1):
        for y in range(x+1,species_num+1,1):
            blast_file=out_folder+'/mSynF'+str(x)+'/species'+str(x)+'_species'+str(y)+'.blast'
            if os.path.exists(blast_file):
                continue
            else:
                ref_prot=out_folder+'/mSynF'+str(x)+'/species.gff_pep_DelTandem'
                query_prot=out_folder+'/mSynF'+str(y)+'/species.gff_pep_DelTandem'
                pick_tandem.do_blast(query_prot,ref_prot,out_folder+'/mSynF'+str(x),tools,threads_num,evalue,'/species'+str(x)+'_species'+str(y)+'.blast')
def get_gene_pos_dict(ref_pos):
    return_list_dict={}
    return_dict={}
    op_pos=open(ref_pos,'r')
    # pos_index=0
    for line in op_pos:
        line_list=line.strip().split('\t')
        # print(line_list)
        if line_list[3] in return_list_dict.keys():
            pos_index+=1
            #[pos_index,length,chr_num]
            return_dict[line_list[0]]=[pos_index,int(line_list[2])/3,line_list[3]]
            return_list_dict[line_list[3]].append(line_list[0])
        else:
            pos_index=60
            # index length chromosome
            return_dict[line_list[0]]=[pos_index,int(line_list[2])/3,line_list[3]]
            # 给加60个空格避免每个染色体前基因不够遍历
            return_list_dict[line_list[3]]=['']*60
            return_list_dict[line_list[3]].append(line_list[0])
    for key,value in return_list_dict.items():
        # 给加60个空格避免每个染色体前基因不够遍历
        return_list_dict[key]=value+['']*60
    op_pos.close()
    return return_dict,return_list_dict
# def flanking(blast_dict):

def synteny(blast_file,ref_pos,query_pos,out_folder):
    op_blast_file=open(blast_file,'r')
    query_pos_dict,query_list_dict=get_gene_pos_dict(query_pos)
    ref_pos_dict,ref_list_dict=get_gene_pos_dict(ref_pos)
    ref_blast_dict={}
    ref_blast_best_dict={}
    query_blast_dict_tmp={}
    query_blast_best_dict={}
    for line in op_blast_file:
        line_list=line.strip().split('\t')
        if line_list[0] in ref_blast_dict.keys():
            ############ 最多只用15同源基因 #############
            if len( ref_blast_dict[line_list[0]])>15:
                continue
            ref_blast_dict[line_list[0]].append([line_list[1],float(line_list[2]),float(line_list[3]),float(line_list[-2])])
        else:
            ref_blast_best_dict[line_list[0]]=line_list[1]
            ref_blast_dict[line_list[0]]=[[line_list[1],float(line_list[2]),float(line_list[3]),float(line_list[-2])]]
        ##### get query best gene and query homologous dict###
        if line_list[1] in query_blast_dict_tmp.keys():
            query_blast_dict_tmp[line_list[1]].append([line_list[0],float(line_list[2]),float(line_list[3]),float(line_list[-2])])
        else:
            query_blast_dict_tmp[line_list[1]]=[[line_list[0],float(line_list[2]),float(line_list[3]),float(line_list[-2])]]
    op_blast_file.close()
    #对query 排序
    query_blast_dict={}
    for key,value in query_blast_dict_tmp.items():
        query_blast_best_dict[key]=sorted(value,key=lambda x: x[1])[0][0]
        # 最多只用15同源基因
        query_blast_dict[key]=sorted(value,key=lambda x: x[1])[:15]
    # loop blast dict
    def get_syn_dict(blast_dict,pos1_dict,pos2_dict,sp1_list_dict,sp2_list_dict,best_blast_dict):
        ref_syn_dict={}
        for key,value in blast_dict.items():
            ref_index=pos1_dict[key][0]
            ref_chr=pos1_dict[key][2]
            ref_length=pos1_dict[key][1]
            key_list=sp1_list_dict[ref_chr][ref_index-60:ref_index]+sp1_list_dict[ref_chr][ref_index+1:ref_index+61]
            key_best_list=[]
            for key_gene in key_list:
                if key_gene=='':
                    continue
                try:
                    key_best_list.append(best_blast_dict[key_gene] )
                except KeyError:
                    continue
            # key_best_list=[ref_blast_best_dict[key_gene] for key_gene in key_list]
            #[pos_index,length,chr_num]
            for query_gene_list in value:
                query_length=pos2_dict[query_gene_list[0]][1]
                ######### evalue #######
                if query_gene_list[3]> 1e-5:
                    continue
                ######### identity #######
                elif query_gene_list[1]<30:
                    continue
                ######### coverage #######
                elif query_gene_list[2]/query_length < 0.3 or query_gene_list[2]/ref_length < 0.3:
                    continue
                query_gene=query_gene_list[0]
                query_chr=pos2_dict[query_gene][2]
                query_index=pos2_dict[query_gene][0]
                up_query_list=sp2_list_dict[query_chr][query_index-20:query_index]
                down_query_list=sp2_list_dict[query_chr][query_index+1:query_index+21]
                up_num=len(list(set(key_best_list) & set(up_query_list)))
                down_num=len(list(set(key_best_list) & set(down_query_list)))
                if up_num > 3 or down_num > 3:
                    if key in ref_syn_dict.keys():
                        ref_syn_dict[key].append([query_gene_list[0],query_gene_list[3],up_num,down_num])
                    else:
                        ref_syn_dict[key]=[[query_gene_list[0],query_gene_list[3],up_num,down_num]]
        return ref_syn_dict
    # check which synfile have more synteic gene number
    def get_dict_num(count_dict):
        total_num=0
        for key,value in count_dict.items():
            total_num+=len(value)
        return total_num

    ref_syn_dict1=get_syn_dict(ref_blast_dict,ref_pos_dict,query_pos_dict,ref_list_dict,query_list_dict,ref_blast_best_dict)
    syn_num1=get_dict_num(ref_syn_dict1)
    ref_syn_dict2=get_syn_dict(query_blast_dict,query_pos_dict,ref_pos_dict,query_list_dict,ref_list_dict,query_blast_best_dict)
    syn_num2=get_dict_num(ref_syn_dict2)
    if syn_num1>syn_num2:
        ref_syn_dict=ref_syn_dict1
    else:
        ref_syn_dict=ref_syn_dict2
    out_file=open(out_folder,'w')
    for key,value in ref_syn_dict.items():
        for i_list in value:
            out_file.write(key+'\t')
            for i in i_list:
                out_file.write(str(i)+'\t')
            out_file.write('\n')   
    out_file.close() 

def run_synteny(out_fold,threads_num):
    pool = multiprocessing.Pool(threads_num)
    for mSynF in os.listdir(out_fold):
        if mSynF.startswith('mSynF'):
            for blast_file in os.listdir(out_fold+'/'+mSynF):
                if blast_file.endswith('.blast'):
                    species1_fold=out_fold+'/mSynF'+blast_file[:-6].split('_')[0][7:]
                    species2_fold=out_fold+'/mSynF'+blast_file[:-6].split('_')[1][7:]
                    loop_blast=out_fold+'/'+mSynF+'/'+blast_file
                    out_syn_file=out_fold+'/'+mSynF+'/'+blast_file[:-6]+'.syn'
                    # print(species1_fold,species2_fold,loop_blast,out_syn_file)
                    result=pool.apply_async(synteny,(loop_blast,species2_fold+'/species.gff_pos_DelTandem',species1_fold+'/species.gff_pos_DelTandem',out_syn_file))
                    # synteny(loop_blast,species2_fold+'/species.gff_pos_DelTandem',species1_fold+'/species.gff_pos_DelTandem',out_syn_file)
                    # syntenic_segment.syntenic_segment_func(species2_fold+'/species.gff_pos_DelTandem',species1_fold+'/species.gff_pos_DelTandem',out_syn_file,out_syn_file+'.seg')
    result.get()
    pool.close()
    pool.join()

    pool = multiprocessing.Pool(threads_num)
    for mSynF in os.listdir(out_fold):
        if mSynF.startswith('mSynF'):
            for blast_file in os.listdir(out_fold+'/'+mSynF):
                if blast_file.endswith('.blast'):
                    species1_fold=out_fold+'/mSynF'+blast_file[:-6].split('_')[0][7:]
                    species2_fold=out_fold+'/mSynF'+blast_file[:-6].split('_')[1][7:]
                    out_syn_file=out_fold+'/'+mSynF+'/'+blast_file[:-6]+'.syn'
                    result=pool.apply_async(syntenic_segment.syntenic_segment_func,(species2_fold+'/species.gff_pos_DelTandem',species1_fold+'/species.gff_pos_DelTandem',out_syn_file,out_syn_file+'.seg'))
                    # syntenic_segment.syntenic_segment_func(species2_fold+'/species.gff_pos_DelTandem',species1_fold+'/species.gff_pos_DelTandem',out_syn_file,out_syn_file+'.seg')
    result.get()
    pool.close()
    pool.join()



if __name__ == "__main__":
    import sys
    # together_blast(sys.argv[1],'diamond',20,'1e-5')
    # synteny(blast_file,ref_pos,query_pos,out_folder)
    # synteny(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    run_synteny(sys.argv[1])

