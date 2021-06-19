# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 18:50:33 2020

@author: OOO
"""
'''

'''
import pandas as pd
import numpy as np


def jc(U1,U2):

    a=list(set(U1).intersection(set(U2)))


    b=list(set(U1).union(set(U2)))

    result=len(a)/len(b)
    return result
    

#circRNA_name
circRNA= pd.read_csv('circRNA_Name.txt', header=None)
circRNA =circRNA.drop_duplicates()
circRNA = np.array(circRNA)
c_Name = [item for item in list(circRNA[:,0])]


fcirc2dis = open('.\circRNA-miRNA\circ2mir_copy.txt')

circ2dis = pd.read_csv(fcirc2dis, sep='\t',header=None)
circ2dis.columns=['circRNA','miRNA']
circRNA =circ2dis['circRNA'].values.tolist()
miRNA =circ2dis['miRNA'].values.tolist()


circ_ID = [x for x in c_Name if x in circRNA] 

circRNA_m = [[] for i in range(len(circ_ID))] 

for i in range(len(circRNA_m)):
    for j in range(len(circRNA)):
        if circ_ID[i]==circRNA[j]:
            circRNA_m[i].append(miRNA[j])

#267*267
circ_sim=[[[] for i in range(len(circRNA_m))] for j in range(len(circRNA_m))]

#JC Score
for i in range(len(circRNA_m)):
    #circ_sim1.append(ssmodel.wv.similarity(circRNA_ID[0],circRNA_ID[i]))
    for j in range(len(circRNA_m)):
        circ_sim[i][j].append(jc(circRNA_m[i],circRNA_m[j]))
test=pd.DataFrame(columns=circ_ID,index=circ_ID,data=circ_sim)
test.to_csv('e:/varify.csv',encoding='gbk')