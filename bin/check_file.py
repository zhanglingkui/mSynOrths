# -*- coding=utf-8 -*-
# 2020.1007
# @zlk
# 检测文件是否存在，并初始化工作文件目录
import sys
import os
def check_file(fa_gff_file,out_folder,file_index):
    species_folder=open(out_folder+'/Species_folder.txt','a+')
    for line in open(fa_gff_file, 'r'):
        line_list = line.strip().split(',')
        if line.strip()=='':
            continue
        elif len(line_list) == 2:
            file_index+=1
            os.mkdir(out_folder+'/mSynF'+str(file_index))
            if os.path.exists(line_list[0]) and os.path.exists(line_list[1]):
                species_folder.write('mSynF'+str(file_index)+'\t'+line)
                fa_abs = os.path.abspath(line_list[0])
                os.system('ln -s '+fa_abs+' '+ out_folder+'/mSynF'+str(file_index)+'/species.fa')
                gff_abs = os.path.abspath(line_list[1])
                os.system('ln -s '+gff_abs+' '+ out_folder+'/mSynF'+str(file_index)+'/species.gff')
            else:
                print(line+' file not exit')
                sys.exit()
        else:
            print('File location file format error ')
            sys.exit()
if __name__ == "__main__":
    check_file(sys.argv[1],sys.argv[2],0)
        