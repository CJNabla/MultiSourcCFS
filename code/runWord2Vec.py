import gensim
import gensim.models
import os
import numpy as np
class MySentences(object):
    def __init__(self, dirname):
        self.dirname = dirname

    def __iter__(self):
        #for line in open('AllAxioms1.lst'):
        for line in open('AllAxioms.lst'):
            yield line.split()

#sentences =gensim.models.word2vec.LineSentence('AllAxioms1.lst') # a memory-friendly iterator

sentences =gensim.models.word2vec.LineSentence('AllAxioms.lst') # a memory-friendly iterator

ssmodel =gensim.models.Word2Vec(sentences,min_count=0, size=200, window=10, sg=1, negative=4, iter=5)
#Store vector of each  class
GOvectors={}
word_vectors=ssmodel.wv
file= open ('VecResults.lst', 'w')
with open('AllClasses.lst') as f:
    for line in f:
        GO_class=line.rstrip()
        if GO_class in word_vectors.vocab:
            GOvectors[GO_class]=ssmodel[GO_class];
            file.write (str(GO_class) + ' '+ str(GOvectors[GO_class]) +'\n')
file.close()



import pandas as pd


# 打开文件
#读取gene_symbol
gsName= pd.read_csv('genesymbol_Name.txt', header=None)
gsName.columns=['gene_symbol']

gsName = gsName['gene_symbol'].drop_duplicates()
gsName = np.array(gsName)
gsName = [item for item in list(gsName)]


circ_sim=[[[] for i in range(len(gsName))] for j in range(len(gsName))]

for i in range(len(gsName)):

    for j in range(len(gsName)):
        circ_sim[i][j].append(ssmodel.wv.similarity(gsName[i],gsName[j]))
test1=pd.DataFrame(columns=gsName,index=gsName,data=circ_sim)

test.to_csv('GeneSymbol_remove.csv',encoding='gbk')