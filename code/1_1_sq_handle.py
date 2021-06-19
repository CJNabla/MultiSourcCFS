# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 07:43:04 2020

@author: OOO
"""
import numpy as np
import pandas as pd
from math import sqrt
#显示所有列
pd.set_option('display.max_columns', None)

def readfile():
    fr=open('sequence.fasta', 'r')

    seq=[]
    for line in fr:
      seq.append(line)
    fr.close()     

    rna_seq=seq
    return rna_seq

def pigment(rna_seq):

    x=[]
    y=[]
    pre_x=0.5
    pre_y=0.5

    size = 2 ** 3
    for i in rna_seq:

        if i=='A':
            pre_x = pre_x * 0.5
            pre_y = pre_y*0.5
            x.append(pre_x)
            y.append(pre_y)
            continue

        if i=='T':
            pre_x = pre_x * 0.5+0.5
            pre_y = pre_y*0.5
            x.append(pre_x)
            y.append(pre_y)
            continue

        if i == 'C':
            pre_x = pre_x * 0.5
            pre_y = pre_y * 0.5 +0.5
            x.append(pre_x)
            y.append(pre_y)
            continue

        if i == 'G':
            pre_x = pre_x * 0.5 +0.5
            pre_y = pre_y * 0.5 +0.5
            x.append(pre_x)
            y.append(pre_y)
            continue
#    print(x)
#    a_x=[]
#    a_y = []
    vector_x=[]
    vector_y = []
    vector_sum=[]
    label=[]

#    print('x='+str(x))

    for i in range(len(x)):

        label.append(int(x[i] * size)+int(y[i] * size)*size)
#    print(label)
    for i in range (size*size):
        x_sum = 0
        y_sum = 0
        sum = 0
        for j in range(len(label)):

            if i == label[j]:

                x_sum = x_sum+x[j]
                y_sum = y_sum +y[j]

                sum += 1


        if sum >0:
            x_sum = x_sum
            y_sum = y_sum
        else:
            x_sum = 0
            y_sum = 0

#        a_x.append(x_sum)
#        a_y.append(y_sum)
        vector_x.append(x_sum)
        vector_y.append(y_sum)
        vector_sum.append(sum)
    

    vector=[]
#    print('sum=\n')
#   print(sum)
#    print('vector_sum=\n')
#   print(len(vector_sum)  )  
    for i in vector_sum:

        if i==0:#sum or i
            vector.append(0)

        else:
            z_score = (i-np.average(vector_sum))/np.std(vector_sum)

            vector.append(z_score)

    for i in vector_x:
        vector.append(i)

    for i in vector_y:
        vector.append(i)

    
    return vector

def P(a,b):
    s1 = pd.Series(a)
    s2 = pd.Series(b)
    corr = s1.corr(s2)

    return corr
def Global_distance(a,b):

    nw = 0
    xw = 0
    yw = 0
    sx = 0
    sy = 0
    rw = 0
    for i in range(len(a)):



        nw = nw+a[i]*b[i]

        xw = xw+a[i]*a[i]*b[i]
        yw = yw+a[i]*b[i]*b[i]
    if nw == 0:
        print(a)
        print(b)

    xw = xw/nw

    yw = yw/nw


    for i in range(len(a)):
        sx=sx+(a[i]-xw)*(a[i]-xw)*a[i]*b[i]
        sy=sy+(b[i]-yw)*(b[i]-yw)*a[i]*b[i]
    if sx==0:

        print(a)
        print(b)
    if sy == 0:
        print("error2")
    sx = sx/nw
    sy = sy/nw



    for i in range(len(a)):

        rw = rw+((a[i]-xw)/sqrt(sx)) * ((b[i]-yw)/sqrt(sy))*a[i]*b[i]

    rw =1 - rw/nw


    return rw
if __name__ == "__main__":


 #   rna_seq=['AC']
    rna_seq=readfile()
    j=0

    rna_name=[]
    seq_data=[]


    rw=[]

    for i in rna_seq:


        if i[1]!='':
           j=j+1

        seq_data.append(pigment(i))

#    print('j',j)

    k=0


    rw1=[[[] for i in range(365)] for j in range(365)]

    for i in range(len(seq_data)):

        rw_part = []
        for j in range(len(seq_data)):

                #rw_part.append( P(i, j) )
                #rw_part.append( Global_distance(i, j) )
                if i!=j:
                    rw1[i][j].append(P(seq_data[i],seq_data[j]))
                else:
                    rw1[i][j].append(1)
                k=k+1

        #rw.append(rw_part)
    
test1=pd.DataFrame(data=rw1)
test1.to_csv('e:/chaos_game_sequece.csv',encoding='gbk')

#    print(k)
