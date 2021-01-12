#!/home/pyenv/versions/py3.7/bin/python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import bisect
import sys
import os

predictor_code = sys.argv[1]
input_path = os.path.abspath(sys.argv[2])
output_path = os.path.abspath(sys.argv[3])

# In[2]:


Lines= []
f = open(input_path,'r')
for line in f:
    data_line = line.rstrip().split('\t')
    Lines.append(data_line)


# In[3]:


Uniprot_ID=[]
AA_Sequence=[]
PredictionScore_byProt=[]
QualityScore_byProt=[]


J = len(Lines)-3
b =0
for b in range(0,J,4):
    Uniprot_ID.append(Lines[b][0])
    AA_Sequence.append(Lines[b+1][0])
    PredictionScore_byProt.append(list(map(float,Lines[b+2][0].rstrip().split(','))))
    QualityScore_byProt.append(list(map(float,Lines[b+3][0].rstrip().split(','))))



# In[4]:


#Getting length of each protein to a list
Length_EachProtein = []
b=0
for b in range(0,len(AA_Sequence), 1):
                   c = len(AA_Sequence[b])
                   Length_EachProtein.append(c)

#Getting continous sum off 999 proteins to an array
d=np.cumsum(Length_EachProtein)


# In[5]:


PredictionScore=np.concatenate(PredictionScore_byProt)
QualityScore=np.concatenate(QualityScore_byProt)




# In[6]:


file = (r"MCC_Lookup.csv")
MCC_Lookup= pd.read_csv(file)
MCC_Lookup.head()


# In[7]:


dis465_QAWindow_tr=MCC_Lookup.iloc[:,(1)].values
dis465_MCCWindow_tr=MCC_Lookup.iloc[:,(2)].values
disHL_QAWindow_tr=MCC_Lookup.iloc[:,(3)].values
disHL_MCCWindow_tr=MCC_Lookup.iloc[:,(4)].values
espD_QAWindow_tr=MCC_Lookup.iloc[:,(5)].values
espD_MCCWindow_tr=MCC_Lookup.iloc[:,(6)].values
espN_QAWindow_tr=MCC_Lookup.iloc[:,(7)].values
espN_MCCWindow_tr=MCC_Lookup.iloc[:,(8)].values
espX_QAWindow_tr=MCC_Lookup.iloc[:,(9)].values
espX_MCCWindow_tr=MCC_Lookup.iloc[:,(10)].values
glo_QAWindow_tr=MCC_Lookup.iloc[:,(11)].values
glo_MCCWindow_tr=MCC_Lookup.iloc[:,(12)].values
iupl_QAWindow_tr=MCC_Lookup.iloc[:,(13)].values
iupl_MCCWindow_tr=MCC_Lookup.iloc[:,(14)].values
iups_QAWindow_tr=MCC_Lookup.iloc[:,(15)].values
iups_MCCWindow_tr=MCC_Lookup.iloc[:,(16)].values
jronn_QAWindow_tr=MCC_Lookup.iloc[:,(17)].values
jronn_MCCWindow_tr=MCC_Lookup.iloc[:,(18)].values
vsl_QAWindow_tr=MCC_Lookup.iloc[:,(19)].values
vsl_MCCWindow_tr=MCC_Lookup.iloc[:,(20)].values
spot_QAWindow_tr=MCC_Lookup.iloc[:,(21)].values
spot_MCCWindow_tr=MCC_Lookup.iloc[:,(22)].values
disopred_QAWindow_tr=MCC_Lookup.iloc[:,(23)].values
disopred_MCCWindow_tr=MCC_Lookup.iloc[:,(24)].values


# In[8]:


import bisect
mcc_Predicted=[]

b=0
for b in range(0,len(QualityScore),1):
         if predictor_code=='0':
                 try:
                     mcc=dis465_MCCWindow_tr[bisect.bisect_left(dis465_QAWindow_tr, QualityScore[b])]
                 except:
                    mcc=dis465_MCCWindow_tr[(bisect.bisect_left(dis465_QAWindow_tr, QualityScore[b]))-1]
         if predictor_code=='1':
                 try:
                   mcc=disHL_MCCWindow_tr[bisect.bisect_left(disHL_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=disHL_MCCWindow_tr[(bisect.bisect_left(disHL_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='2':		   
                 try:
                   mcc=espD_MCCWindow_tr[bisect.bisect_left(espD_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=espD_MCCWindow_tr[(bisect.bisect_left(espD_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='3':		   
                 try:
                   mcc=espN_MCCWindow_tr[bisect.bisect_left(espN_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=espN_MCCWindow_tr[(bisect.bisect_left(espN_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='4':		   
                 try:
                   mcc=espX_MCCWindow_tr[bisect.bisect_left(espX_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=espX_MCCWindow_tr[(bisect.bisect_left(espX_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='5':		   
                 try:
                   mcc=glo_MCCWindow_tr[bisect.bisect_left(glo_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=glo_MCCWindow_tr[(bisect.bisect_left(glo_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='6':		   
                 try:
                   mcc=iupl_MCCWindow_tr[bisect.bisect_left(iupl_QAWindow_tr, QualityScore[b])]
                 except:
                    mcc=iupl_MCCWindow_tr[(bisect.bisect_left(iupl_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='7':			
                 try:
                   mcc=iups_MCCWindow_tr[bisect.bisect_left(iups_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=iups_MCCWindow_tr[(bisect.bisect_left(iups_QAWindow_tr, QualityScore[b]))-1]
         
         if predictor_code=='9':		   
                 try:
                   mcc=vsl_MCCWindow_tr[bisect.bisect_left(vsl_QAWindow_tr, QualityScore[b])]
                 except:
                   mcc=vsl_MCCWindow_tr[(bisect.bisect_left(vsl_QAWindow_tr, QualityScore[b]))-1]
         
         mcc_Predicted.append(mcc)
       
           


# In[9]:


def Score_Binary(Score,thresh):
                     Binary=[]
                     b=0
                     for b in range(0,len(Score),1):
                            if Score[b]>=thresh:Binary.append(1)
                            else:Binary.append(0)
                     return Binary   


# In[10]:




if predictor_code=='0':
           Binary=Score_Binary(PredictionScore,0.5)
if predictor_code=='1':
           Binary=Score_Binary(PredictionScore,0.086)
if predictor_code=='2':
           Binary=Score_Binary(PredictionScore,0.507)
if predictor_code=='3':
           Binary=Score_Binary(PredictionScore,0.309)
if predictor_code=='4':
           Binary=Score_Binary(PredictionScore,0.143)
if predictor_code=='5':
           Binary=Score_Binary(PredictionScore,0.0)
if predictor_code=='6':
           Binary=Score_Binary(PredictionScore,0.5)
if predictor_code=='7':
           Binary=Score_Binary(PredictionScore,0.5)
if predictor_code=='9':
           Binary=Score_Binary(PredictionScore,0.5)


# In[11]:



Binary_byProt= np.split(Binary,d)




# In[12]:





# In[32]:


MCC_byProt=np.split(mcc_Predicted,d)




# In[33]:


Lines=[]
b=0
for b in range(0,len(Uniprot_ID),1):
                Lines.append(Uniprot_ID[b])
                Lines.append(AA_Sequence[b])
                x = str(list(Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
           
                



# In[34]:


with open(output_path, 'w') as f:
    for item in Lines:
        f.write("%s\n" % item)

