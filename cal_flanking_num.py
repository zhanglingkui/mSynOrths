import sys

file1=open(sys.argv[1],'r')
num_dict={}
num_list=[]
syn_dict={}
for line in file1:
    line_list=line.strip().split('\t')
    if int(line_list[-1])>int(line_list[-2]):
        num_list.append(int(line_list[-1]))
    else:
        num_list.append(int(line_list[-2]))
import matplotlib.pyplot as plt
import seaborn as sns
plt.hist(num_list)
# sns.distplot(num_list,norm_hist=False,kde=False)
plt.savefig('xx.png')