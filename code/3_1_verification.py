# -*- coding: utf-8 -*-
"""

@author: OOO
"""

import pandas as pd
import numpy as np

import xlrd

def import_excel_matrix(path):
  table = xlrd.open_workbook(path).sheets()[0] 
  row = table.nrows 
  col = table.ncols
  datamatrix = np.zeros((row, col)) 
  for i in range(col): 
    cols = np.matrix(table.col_values(i)) 
    datamatrix[:, i] = cols 
  return datamatrix


data_file = u'result.xlsx' 

gs=import_excel_matrix(data_file)

#circRNA_name
circRNA= pd.read_csv('circRNA_Name.txt', header=None)
circRNA =circRNA.drop_duplicates()
circRNA = np.array(circRNA)
c_Name = [item for item in list(circRNA[:,0])]

data = 'circid_m.txt' 
circ_m= pd.read_csv(open(data), header=None)
circ_m =circ_m.drop_duplicates()
circ_m = np.array(circ_m)
c_m = [item for item in list(circ_m[:,0])]



circRNA_m = []

for i in range(len(c_m)):
    for j in range(len(c_Name)):
        if c_m[i]==c_Name[j]:
            circRNA_m.append(j)

#267*267
circ_sim=[[[] for i in range(len(c_m))] for j in range(len(c_m))]


for i in range(len(c_m)):
    #circ_sim1.append(ssmodel.wv.similarity(circRNA_ID[0],circRNA_ID[i]))
    for j in range(len(c_m)):
        circ_sim[i][j].append(gs[circRNA_m[i]][circRNA_m[j]])

test=pd.DataFrame(columns=c_m,index=c_m,data=circ_sim)
test.to_csv('sv.csv',encoding='gbk')

