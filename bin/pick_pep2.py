# -*- coding=utf-8 -*-
# 2021.0708
# @zlk 
# zhanglk960127@163.com
#提取坐标文件和蛋白文件
import os
import sys
def pick_mRNA(gff_file,fa_file,out_folder,gffread_soft):
    
    def get_pep(gffread,out_file):
        os.system("%s %s/species.gff -g %s/species.fa -S -y %s/species.gff_pep" %(gffread,out_file,out_file,out_file))
    get_pep(gffread_soft,out_folder)
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
    def get_pos_file(in_gff_file,out_file):
        op_gff_file=open(in_gff_file,'r')
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
        
        cds_parent_dict={}
        orginal_chr_list=[]
        for i_list in gff_list:
            if i_list[2]!='CDS':
                continue
            if i_list[0] not in orginal_chr_list:
                orginal_chr_list.append(i_list[0])
            parent_str=pick_parent(i_list[8])
            if parent_str in cds_parent_dict.keys():
                cds_parent_dict[parent_str].append([int(i_list[3]),int(i_list[4])])
            else:
                #####  + chr cds_start cds_end #######
                cds_parent_dict[parent_str]=[i_list[6],i_list[0],[int(i_list[3]),int(i_list[4])]]
        for key,value in cds_parent_dict.items():
            total_length=0
            for i in value[3:]:
                total_length+=abs(i[1]-i[0])
            #####  + length chr cds_start cds_end #######
            cds_parent_dict[key].insert(1,total_length)
        mRNA_pos_dict={}
        if gff_type=='gene' or gff_type=='mRNA':
            for i_list in gff_list:
                if i_list[2]=='gene' or i_list=='mRNA':
                    mrna_id=pick_id[i_list[8]]
                    if i_list[0] in mRNA_pos_dict.keys():
                        mRNA_pos_dict[i_list[0]].append([mrna_id,int(i_list[3]),int(i_list[4]),i_list[6]])
                    else:
                        mRNA_pos_dict[i_list[0]][[mrna_id,int(i_list[3]),int(i_list[4]),i_list[6]]]

        else:
            mrna_parent_dict={}
            for i_list in gff_list:
                if i_list[2]!='mRNA':
                    continue
                parent_str=pick_parent(i_list[8])
                mrna_id=pick_id[i_list[8]]
                if i_list[0] in mRNA_pos_dict.keys():
                    mRNA_pos_dict[i_list[0]].append([mrna_id,int(i_list[3]),int(i_list[4]),i_list[6]])
                else:
                    # mRNA_pos_dict[i_list[0]]={}
                    mRNA_pos_dict[i_list[0]][[mrna_id,int(i_list[3]),int(i_list[4]),i_list[6]]]
                # mRNA_pos_dict[mrna_id][i_list[0],int(i_list[3]),int(i_list[4]),i_list[6]]
                if parent_str in mrna_parent_dict.keys():
                    mrna_parent_dict[parent_str].append(mrna_id)
                else:
                    mrna_parent_dict[parent_str]=[mrna_id]
            have_mrna_list=[]
            for key,value in mrna_parent_dict.items():
                if len(value)>1:
                    sort_value=sorted(value,key=lambda x : cds_parent_dict[x][1],reverse=True)
                    have_mrna_list.append(sort_value[0])
                else:
                    have_mrna_list.append(value[0])
            # ###### 修改######
            # for key,value in mRNA_pos_dict.items():
            #     if key not in have_mrna_list:
            #         del mRNA_pos_dict[key]
        out_wr_file=open(out_file,'w')
        for i in orginal_chr_list:
            sort_list=sorted(mRNA_pos_dict[i],key=lambda x : x[1])
            for mrna_list in sort_list:
                if mrna_list[0] in have_mrna_list:
                    out_wr_file.write(mrna_list[0]+'\t'+key+'\t'+str(mrna_list[1])+'\t'+str(mrna_list[2])+'\t'+cds_parent_dict[mrna_list[0]][1]+'\t'+mrna_list[3]+'\n')         

if __name__ == '__main__':
    
    main()
    

        






