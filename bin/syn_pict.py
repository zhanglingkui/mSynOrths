import collections
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype']=42
def get_pos_dict(pos_file,chr_gap):
    loop_pos_file=open(pos_file,'r')
    chr_gene_pos_dict=collections.OrderedDict()
    for line in loop_pos_file:
        line_list=line.strip().split('\t')
        if line_list[3] in chr_gene_pos_dict.keys():
            chr_gene_pos_dict[line_list[3]].append([line_list[0],int(line_list[4].split(':')[0])])
        else:
            chr_gene_pos_dict[line_list[3]]=[[line_list[0],int(line_list[4].split(':')[0])]]

    loop_pos_file.close()
    gene_pos_dict={}
    chr_pos_list=[]
    chr_start=0
    for key,value in chr_gene_pos_dict.items():
        if len(value)<300:
            continue
        # chr_name, chr_start_pos ,chr_end_pos
        for i in value:
            gene_pos_dict[i[0]]=i[1]+chr_start
        chr_pos_list.append([key,chr_start,chr_start+value[-1][1]-value[0][1]])
        chr_start+=(value[-1][1]-value[0][1]+chr_gap)
    return gene_pos_dict,chr_pos_list

def syn_point_pict(pos_file1,pos_file2,syn_file):
    gene_pos_dict1,chr_pos_list1=get_pos_dict(pos_file1,0)
    gene_pos_dict2,chr_pos_list2=get_pos_dict(pos_file2,0)
    fig = plt.figure(figsize=(10, 10))
    scatter_x,scatter_y=[],[]
    loop_syn_file=open(syn_file,'r')
    for line in loop_syn_file:
        if line.startswith('###'):
            continue
        line_list=line.strip().split('\t')
        if line_list[0] in gene_pos_dict1.keys() and line_list[1] in gene_pos_dict2.keys():
            scatter_x.append(gene_pos_dict1[line_list[0]])
            scatter_y.append(gene_pos_dict2[line_list[1]])
    


    loop_syn_file.close()
    plt.scatter(scatter_x,scatter_y,color='red',s=0.1,alpha = 1)
    plt.savefig('xx.png',dpi=100)





if __name__ == '__main__':
    import sys

    syn_point_pict(sys.argv[1],sys.argv[2],sys.argv[3])
    
    


