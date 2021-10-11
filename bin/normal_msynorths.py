# -*- coding=utf-8 -*-
# 2020.1008
# @zlk
# 将正常共线性的流程穿起来
import os
import sys
import multiprocessing
import pick_pep
import check_file
import time
import pick_tandem
import synteny

def normal(fa_gff_file,out_folder,threads,tools):
    # if os.path.exists(out_folder):
    #     if len(os.listdir(out_folder))!=0:
    #         print('There are files in the output folder, please change the output folder!!!')
    #         sys.exit()
    # else:
    #     os.mkdir(out_folder)

    # print(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+'\tStart to initialize the working directory and extract the protein sequence.')
    # check_file.check_file(fa_gff_file,out_folder,0)

    # pick_pep.run_pick_pep(threads,out_folder)

    # print(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+'\tStart to self blast.')
    # pick_tandem.self_blast(out_folder,tools,threads)

    # print(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+'\tStart to pick tandem array.')
    # pick_tandem.run_pick_tandem(out_folder,threads)

    # print(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+'\tStart to together blastp.')
    # synteny.together_blast(out_folder,tools,threads,'1e-5')

    # print(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+'\tStart to syntenic analyse.')
    synteny.run_synteny(out_folder,20)



if __name__ == "__main__":
    normal(sys.argv[1],sys.argv[2],20,sys.argv[3])



