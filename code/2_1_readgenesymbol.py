# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 08:55:26 2020
得到每个genesymbol关联的GO
@author: Solar
"""

import pandas as pd
import numpy as np
#读取数据 
    
#显示所有列
pd.set_option('display.max_columns', None)
#%%
#读取gene_symbol
gsName= pd.read_csv('genesymbol_Name.txt', header=None)
gsName.columns=['gene_symbol']

gsName = gsName['gene_symbol'].drop_duplicates()
gsName = np.array(gsName)
gsName = [item for item in list(gsName)]

#读取symbol2GO数据——10090版本
f_10090 = open('E:\\学习心得\circRNA Functional Similarity\数据\gene_symbol2GO\\_dbOrg_GO_ID__to__Gene_Symbol_10090')

g2g_10090 = pd.read_csv(f_10090, sep='\t',header=None)
g2g_10090.columns=['GO_id','gene_symbol']
GO_10090 =g2g_10090['GO_id'].values.tolist()
GS_10090 =g2g_10090['gene_symbol'].values.tolist()

#g2g_10090=np.array(g2g_10090)
#g2g_10090 = g2g_10090['GO_id'].values.tolist()
#GO_10090 = [item for item in list(g2g_10090[:,0])]
#GS_10090 = [item for item in list(g2g_10090[:,1])]

#读取symbol2GO数据——9606版本
f_9606 = open('E:\\学习心得\circRNA Functional Similarity\数据\gene_symbol2GO\\_dbOrg_GO_ID__to__Gene_Symbol_9606')

g2g_9606 = pd.read_csv(f_9606, sep='\t',header=None)
g2g_9606.columns=['GO_id','gene_symbol']

GO_9606 =g2g_9606['GO_id'].values.tolist()
GS_9606 =g2g_9606['gene_symbol'].values.tolist()
#g2g_10090 = g2g_10090['GO_id'].values.tolist()
#GO_9606 = [item for item in list(g2g_9606[:,0])]
#GS_9606 = [item for item in list(g2g_9606[:,1])]

#读取symbol2GO数据——9544版本
f_9544 = open('E:\\学习心得\circRNA Functional Similarity\数据\gene_symbol2GO\\_dbOrg_GO_ID__to__Gene_Symbol_9544')

g2g_9544 = pd.read_csv(f_9544, sep='\t',header=None)
g2g_9544.columns=['GO_id','gene_symbol']
GO_9544 =g2g_9544['GO_id'].values.tolist()
GS_9544 =g2g_9544['gene_symbol'].values.tolist()
#g2g_9544 =np.array(g2g_9544)
#g2g_10090 = g2g_10090['GO_id'].values.tolist()
#GO_9544 = [item for item in list(g2g_9606[:,0])]
#GS_9544 = [item for item in list(g2g_9606[:,1])]

#读取symbol2GO数据——7955版本
f7955 = open('E:\\学习心得\circRNA Functional Similarity\数据\gene_symbol2GO\\_dbOrg_GO_ID__to__Gene_Symbol_7955')

g2g_7955 = pd.read_csv(f7955, sep='\t',header=None)
g2g_7955.columns=['GO_id','gene_symbol']
GO_7955 =g2g_7955['GO_id'].values.tolist()
GS_7955 =g2g_7955['gene_symbol'].values.tolist()
#g2g_9606 =np.array(g2g_7955)
#g2g_10090 = g2g_10090['GO_id'].values.tolist()
#GO_7955 = [item for item in list(g2g_7955[:,0])]
#GS_7955 = [item for item in list(g2g_7955[:,1])]

#%%得到gene_symbol关联的GO

g2g=[]

G_id1 = [[ 0 for i in range(1)] for row in range(len(gsName))]
G_id2 = [[ 0 for i in range(1)] for row in range(len(gsName))]
G_id3 = [[ 0 for i in range(1)] for row in range(len(gsName))]
G_id4 = [[ 0 for i in range(1)] for row in range(len(gsName))]
for i in range(len(gsName)):
    
    for j in range(len(GO_7955)):
        if gsName[i].lower() in str(GS_7955[j]).lower():
            G_id1[i].append(GO_7955[j])
    for j in range(len(GO_9544)):
        if gsName[i].lower() in GS_9544[j].lower():
             G_id2[i].extend(GO_9544[j])        
    for j in range(len(GO_9606)):
          if gsName[i].lower() in GS_9606[j].lower():
             G_id3[i].append(GO_9606[j])      
    for j in range(len(GO_10090)):
          if gsName[i].lower() in GS_10090[j].lower():
             G_id4[i].append(GO_10090[j])    
G_G=[]

G_G=G_id1+G_id2+G_id3+G_id4
g2g=list(zip(G_id1,G_id2,G_id3,G_id4))

#%%构建语料库

file = open("GO7_7")
fw=open('axioms7_7', 'w')


k=0
for line in file.readlines():
    #print(line)
    name=line.rstrip('\n')
    name=name.split(' ')
    
    for j in range(len(name)):
        fw.write(gsName[k]+' hasFunction '+name[j])
        fw.write('\n')
    k=k+1
    
file.close()


fw.close()
