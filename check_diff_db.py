# -*- coding=utf-8 -*-
# 2021.3.2
# @zlk
# 对不同建库物种选择的结果进行评估

import sys

file1=open(sys.argv[1],'r') #Ath_db
file2=open(sys.argv[2],'r') #Brapa_db
Ath_dict={}
for line in file1:
    line_list=line.strip().split('\t')
    if line_list[0] in Ath_dict.keys():
        Ath_dict[line_list[0]].append(line_list[1])
    else:
        Ath_dict[line_list[0]]=[line_list[1]]
Brapa_dict={}
for line in file2:
    line_list=line.strip().split('\t')
    if line_list[1] in Brapa_dict.keys():
        Brapa_dict[line_list[1]].append([line_list[0],float(line_list[-2]),float(line_list[2])])
    else:
        Brapa_dict[line_list[1]]=[[line_list[0],float(line_list[-2]),float(line_list[2])]]

for key,value in Brapa_dict.items():
    value=sorted(value,key=lambda x : x[2],reverse=False)
    out_list=[x[0] for x in sorted(value,key=lambda x:x[1],reverse=True)]
    
    Brapa_dict[key]=out_list
index,index1,index2=0,0,0
for key,value in Ath_dict.items():
    Ath_value=list(set(value))
    if key in Brapa_dict.keys():
        Brapa_value=list(set(Brapa_dict[key]))
        for i in Ath_value[:9]:
            index1+=1
            if i in Brapa_value[:9]:
                index+=1
        for m in Brapa_value[:9]:
            index2+=1
print(len(Ath_dict))
print(len(Brapa_dict))

print(index,index1,index2)