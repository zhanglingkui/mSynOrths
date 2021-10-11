from sys import path
import synteny
def add_dict(key_str,value_str,new_dict):
    if key_str in new_dict.keys():
        new_dict[key_str].append(value_str)
    else:
        new_dict[key_str]=[value_str]
        
def syntenic_segment_func(pos1,pos2,syn_file,out_seg_file):
    ##_pos_dict:index length chromosome
    syn_dict={}
    syn_list=[]
    homo_dict1={}
    homo_dict2={}
    loop_syn_file=open(syn_file,'r')
    for line in loop_syn_file:
        line_list=line.strip().split('\t')
        add_dict(line_list[0],[line_list[1],float(line_list[2])],homo_dict1)
        add_dict(line_list[1],[line_list[0],float(line_list[2])],homo_dict2)
    loop_syn_file.close()
    loop_syn_file=open(syn_file,'r')
    print(1111)
    print(len(homo_dict1))
    print(len(homo_dict2))
    if len(homo_dict1)>len(homo_dict2):
        homo_dict=homo_dict1
        ref_pos_dict,ref_gene_list=synteny.get_gene_pos_dict(pos1)
        query_pos_dict,query_gene_list=synteny.get_gene_pos_dict(pos2)
        for line in loop_syn_file:
            line_list=line.strip().split('\t')
            add_dict(line_list[0],line_list[1],syn_dict)
            syn_list.append(line_list[:2])
    else:
        homo_dict=homo_dict2
        ref_pos_dict,ref_gene_list=synteny.get_gene_pos_dict(pos2)
        query_pos_dict,query_gene_list=synteny.get_gene_pos_dict(pos1)
        for line in loop_syn_file:
            line_list=line.strip().split('\t')
            add_dict(line_list[1],line_list[0],syn_dict)
            syn_list.append([line_list[1],line_list[0]])
            
    loop_syn_file.close()
    syn_seg_list=[]
    print(1111)
    seg_gap_distance=50
    for syn_pair in syn_list:
        ref_pos=ref_pos_dict[syn_pair[0]]
        query_pos=query_pos_dict[syn_pair[1]]
        add_seg_list=[syn_pair]
        for i in ref_gene_list[ref_pos[2]][ref_pos[0]+1:]:
            if i=='':
                continue
            if i in syn_dict.keys() and abs(ref_pos_dict[i][0]-ref_pos[0])<seg_gap_distance and ref_pos_dict[i][2]==ref_pos[2]:
                for j in syn_dict[i]:
                    if abs(query_pos_dict[j][0]-query_pos[0])<seg_gap_distance and query_pos_dict[j][2]==query_pos[2]:
                        ref_pos=ref_pos_dict[i]
                        query_pos=query_pos_dict[j]
                        if [i,j] in syn_list:
                            syn_list.remove([i,j])
                            add_seg_list.append([i,j])
            elif abs(ref_pos_dict[i][0]-ref_pos[0])>seg_gap_distance:
                break
        ref_pos=ref_pos_dict[syn_pair[0]]
        query_pos=query_pos_dict[syn_pair[1]]
        for i in reversed(ref_gene_list[ref_pos[2]][:ref_pos[0]]):
            if i=='':
                continue
            if i in syn_dict.keys() and abs(ref_pos_dict[i][0]-ref_pos[0])<seg_gap_distance and ref_pos_dict[i][2]==ref_pos[2]:
                for j in syn_dict[i]:
                    if abs(query_pos_dict[j][0]-query_pos[0])<seg_gap_distance and query_pos_dict[j][2]==query_pos[2]:
                        ref_pos=ref_pos_dict[i]
                        query_pos=query_pos_dict[j]
                        if [i,j] in syn_list:
                            syn_list.remove([i,j])
                            add_seg_list.insert(0,[i,j])
            elif abs(ref_pos_dict[i][0]-ref_pos[0])>seg_gap_distance:
                break
        if len(add_seg_list)>10:
            syn_seg_list.append(add_seg_list)
    # print(syn_seg_list)
    out_file=open(out_seg_file,'w')
    for seg_list in syn_seg_list:
        best_hit_num=0
        for gene_pair in seg_list:
            sort_homo_list=sorted(homo_dict[gene_pair[0]],key=lambda x:x[1])
            if gene_pair[1]==sort_homo_list[0][0]:
                best_hit_num+=1
        ### best-hit gene ration #####
        if best_hit_num/len(seg_list)<0.3:
            continue
        out_file.write('######'+str(best_hit_num)+'\n')
    
        for gene_pair in seg_list:
            out_file.write(gene_pair[0]+'\t'+gene_pair[1]+'\n')
    #########根据同源besthit基因数量决定并去除不属于他的共线性片段#####   


if __name__ == '__main__':
    import sys
    syntenic_segment(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

    



        
    




