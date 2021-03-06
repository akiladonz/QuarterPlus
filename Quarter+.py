#!/home/pyenv/versions/py3.7/bin/python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import bisect
import keras
import sys
import tensorflow as tf
from keras import metrics


# In[2]:


Lines= []
f = open(r'Qplus_Input.txt')
for line in f:
    data_line = line.rstrip().split('\t')
    Lines.append(data_line)


# In[3]:


Uniprot_ID=[]
AA_Sequence=[]
dis465_PredictionScore_byProt=[]
dis465_QualityScore_byProt=[]
disHL_PredictionScore_byProt=[]
disHL_QualityScore_byProt=[]
espD_PredictionScore_byProt=[]
espD_QualityScore_byProt=[]
espN_PredictionScore_byProt=[]
espN_QualityScore_byProt=[]
espX_PredictionScore_byProt=[]
espX_QualityScore_byProt=[]
glo_PredictionScore_byProt=[]
glo_QualityScore_byProt=[]
iupl_PredictionScore_byProt=[] 
iupl_QualityScore_byProt=[]
iups_PredictionScore_byProt=[]
iups_QualityScore_byProt=[]
vsl_PredictionScore_byProt=[]
vsl_QualityScore_byProt=[]
spot_PredictionScore_byProt=[]
spot_QualityScore_byProt=[]
disopred_PredictionScore_byProt=[]
disopred_QualityScore_byProt=[]
J = len(Lines)-23
b =0
for b in range(0,J,24):
    Uniprot_ID.append(Lines[b][0])
    AA_Sequence.append(Lines[b+1][0])
    dis465_PredictionScore_byProt.append(list(map(float,Lines[b+2][0].rstrip().split(','))))
    dis465_QualityScore_byProt.append(list(map(float,Lines[b+3][0].rstrip().split(','))))
    disHL_PredictionScore_byProt.append(list(map(float,Lines[b+4][0].rstrip().split(','))))
    disHL_QualityScore_byProt.append(list(map(float,Lines[b+5][0].rstrip().split(','))))
    espD_PredictionScore_byProt.append(list(map(float,Lines[b+6][0].rstrip().split(','))))
    espD_QualityScore_byProt.append(list(map(float,Lines[b+7][0].rstrip().split(','))))
    espN_PredictionScore_byProt.append(list(map(float,Lines[b+8][0].rstrip().split(','))))
    espN_QualityScore_byProt.append(list(map(float,Lines[b+9][0].rstrip().split(','))))
    espX_PredictionScore_byProt.append(list(map(float,Lines[b+10][0].rstrip().split(','))))
    espX_QualityScore_byProt.append(list(map(float,Lines[b+11][0].rstrip().split(','))))
    glo_PredictionScore_byProt.append(list(map(float,Lines[b+12][0].rstrip().split(','))))
    glo_QualityScore_byProt.append(list(map(float,Lines[b+13][0].rstrip().split(','))))
    iupl_PredictionScore_byProt.append(list(map(float,Lines[b+14][0].rstrip().split(','))))
    iupl_QualityScore_byProt.append(list(map(float,Lines[b+15][0].rstrip().split(','))))
    iups_PredictionScore_byProt.append(list(map(float,Lines[b+16][0].rstrip().split(','))))
    iups_QualityScore_byProt.append(list(map(float,Lines[b+17][0].rstrip().split(','))))
    vsl_PredictionScore_byProt.append(list(map(float,Lines[b+18][0].rstrip().split(','))))
    vsl_QualityScore_byProt.append(list(map(float,Lines[b+19][0].rstrip().split(','))))
    spot_PredictionScore_byProt.append(list(map(float,Lines[b+20][0].rstrip().split(','))))
    spot_QualityScore_byProt.append(list(map(float,Lines[b+21][0].rstrip().split(','))))
    disopred_PredictionScore_byProt.append(list(map(float,Lines[b+22][0].rstrip().split(','))))
    disopred_QualityScore_byProt.append(list(map(float,Lines[b+23][0].rstrip().split(','))))


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


dis465_PredictionScore=np.concatenate(dis465_PredictionScore_byProt)
dis465_QualityScore=np.concatenate(dis465_QualityScore_byProt)
disHL_PredictionScore=np.concatenate(disHL_PredictionScore_byProt)
disHL_QualityScore=np.concatenate(disHL_QualityScore_byProt)
espD_PredictionScore=np.concatenate(espD_PredictionScore_byProt)
espD_QualityScore=np.concatenate(espD_QualityScore_byProt)
espN_PredictionScore=np.concatenate(espN_PredictionScore_byProt)
espN_QualityScore=np.concatenate(espN_QualityScore_byProt)
espX_PredictionScore=np.concatenate(espX_PredictionScore_byProt)
espX_QualityScore=np.concatenate(espX_QualityScore_byProt)
glo_PredictionScore=np.concatenate(glo_PredictionScore_byProt)
glo_QualityScore= np.concatenate(glo_QualityScore_byProt)
iupl_PredictionScore=np.concatenate(iupl_PredictionScore_byProt)
iupl_QualityScore=np.concatenate(iupl_QualityScore_byProt)
iups_PredictionScore=np.concatenate(iups_PredictionScore_byProt)
iups_QualityScore=np.concatenate(iups_QualityScore_byProt)
vsl_PredictionScore=np.concatenate(vsl_PredictionScore_byProt)
vsl_QualityScore= np.concatenate(vsl_QualityScore_byProt)
spot_PredictionScore=np.concatenate(spot_PredictionScore_byProt)
spot_QualityScore=np.concatenate(spot_QualityScore_byProt)
disopred_PredictionScore=np.concatenate(disopred_PredictionScore_byProt)
disopred_QualityScore=np.concatenate(disopred_QualityScore_byProt)


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
dis465_mcc_Predicted=[]
disHL_mcc_Predicted=[]
espD_mcc_Predicted=[]
espN_mcc_Predicted=[]
espX_mcc_Predicted=[]
glo_mcc_Predicted=[]
iupl_mcc_Predicted=[]
iups_mcc_Predicted=[]
jronn_mcc_Predicted=[]
vsl_mcc_Predicted=[]
spot_mcc_Predicted=[]
disopred_mcc_Predicted=[]
b=0
for b in range(0,len(dis465_QualityScore),1):
                 try:
                   dis465_mcc=dis465_MCCWindow_tr[bisect.bisect_left(dis465_QAWindow_tr, dis465_QualityScore[b])]
                 except:
                   dis465_mcc=dis465_MCCWindow_tr[(bisect.bisect_left(dis465_QAWindow_tr, dis465_QualityScore[b]))-1]
                 try:
                   disHL_mcc=disHL_MCCWindow_tr[bisect.bisect_left(disHL_QAWindow_tr, disHL_QualityScore[b])]
                 except:
                   disHL_mcc=disHL_MCCWindow_tr[(bisect.bisect_left(disHL_QAWindow_tr, disHL_QualityScore[b]))-1]
                 try:
                   espD_mcc=espD_MCCWindow_tr[bisect.bisect_left(espD_QAWindow_tr, espD_QualityScore[b])]
                 except:
                   espD_mcc=espD_MCCWindow_tr[(bisect.bisect_left(espD_QAWindow_tr, espD_QualityScore[b]))-1]
                 try:
                   espN_mcc=espN_MCCWindow_tr[bisect.bisect_left(espN_QAWindow_tr, espN_QualityScore[b])]
                 except:
                   espN_mcc=espN_MCCWindow_tr[(bisect.bisect_left(espN_QAWindow_tr, espN_QualityScore[b]))-1]
                 try:
                   espX_mcc=espX_MCCWindow_tr[bisect.bisect_left(espX_QAWindow_tr, espX_QualityScore[b])]
                 except:
                   espX_mcc=espX_MCCWindow_tr[(bisect.bisect_left(espX_QAWindow_tr, espX_QualityScore[b]))-1]
                 try:
                   glo_mcc=glo_MCCWindow_tr[bisect.bisect_left(glo_QAWindow_tr, glo_QualityScore[b])]
                 except:
                   glo_mcc=glo_MCCWindow_tr[(bisect.bisect_left(glo_QAWindow_tr, glo_QualityScore[b]))-1]
                 try:
                   iupl_mcc=iupl_MCCWindow_tr[bisect.bisect_left(iupl_QAWindow_tr, iupl_QualityScore[b])]
                 except:
                    iupl_mcc=iupl_MCCWindow_tr[(bisect.bisect_left(iupl_QAWindow_tr, iupl_QualityScore[b]))-1]
                 try:
                   iups_mcc=iups_MCCWindow_tr[bisect.bisect_left(iups_QAWindow_tr, iups_QualityScore[b])]
                 except:
                   iups_mcc=iups_MCCWindow_tr[(bisect.bisect_left(iups_QAWindow_tr, iups_QualityScore[b]))-1]
                 try:
                   vsl_mcc=vsl_MCCWindow_tr[bisect.bisect_left(vsl_QAWindow_tr, vsl_QualityScore[b])]
                 except:
                   vsl_mcc=vsl_MCCWindow_tr[(bisect.bisect_left(vsl_QAWindow_tr, vsl_QualityScore[b]))-1]
                 try:
                   spot_mcc=spot_MCCWindow_tr[bisect.bisect_left(spot_QAWindow_tr, spot_QualityScore[b])]
                 except:
                    spot_mcc=spot_MCCWindow_tr[(bisect.bisect_left(spot_QAWindow_tr, spot_QualityScore[b]))-1]
                 try:
                   disopred_mcc=disopred_MCCWindow_tr[bisect.bisect_left(disopred_QAWindow_tr, disopred_QualityScore[b])]
                 except:
                    disopred_mcc=disopred_MCCWindow_tr[(bisect.bisect_left(disopred_QAWindow_tr, disopred_QualityScore[b]))-1]
                 dis465_mcc_Predicted.append(dis465_mcc)
                 disHL_mcc_Predicted.append(disHL_mcc)
                 espD_mcc_Predicted.append(espD_mcc)
                 espN_mcc_Predicted.append(espN_mcc)
                 espX_mcc_Predicted.append(espX_mcc)
                 glo_mcc_Predicted.append(glo_mcc)
                 iupl_mcc_Predicted.append(iupl_mcc)
                 iups_mcc_Predicted.append(iups_mcc)
                 vsl_mcc_Predicted.append(vsl_mcc)
                 spot_mcc_Predicted.append(spot_mcc)
                 disopred_mcc_Predicted.append(disopred_mcc)            


# In[9]:


def Score_Binary(Score,thresh):
                     Binary=[]
                     b=0
                     for b in range(0,len(Score),1):
                            if Score[b]>=thresh:Binary.append(1)
                            else:Binary.append(0)
                     return Binary   


# In[10]:


disopred_Binary=Score_Binary(disopred_PredictionScore,0.5)
spot_Binary=Score_Binary(spot_PredictionScore,0.49)
iups_Binary=Score_Binary(iups_PredictionScore,0.5)

dis465_Binary=Score_Binary(disopred_PredictionScore,0.5)
disHL_Binary=Score_Binary(spot_PredictionScore,0.086)
espD_Binary=Score_Binary(iups_PredictionScore,0.507)

espN_Binary=Score_Binary(disopred_PredictionScore,0.309)
espX_Binary=Score_Binary(spot_PredictionScore,0.143)
glo_Binary=Score_Binary(iups_PredictionScore,0.0)

iupl_Binary=Score_Binary(disopred_PredictionScore,0.5)
vsl_Binary=Score_Binary(spot_PredictionScore,0.5)


# In[11]:


disopred_Binary_byProt= np.split(disopred_Binary,d)
spot_Binary_byProt= np.split(spot_Binary,d)
iups_Binary_byProt= np.split(iups_Binary,d)

dis465_Binary_byProt=np.split(dis465_Binary,d)
disHL_Binary_byProt=np.split(disHL_Binary,d)
espD_Binary_byProt=np.split(espD_Binary,d)

espN_Binary_byProt=np.split(espN_Binary,d)
espX_Binary_byProt=np.split(espX_Binary,d)
glo_Binary_byProt=np.split(glo_Binary,d)

iupl_Binary_byProt=np.split(iupl_Binary,d)
vsl_Binary_byProt=np.split(vsl_Binary,d)


# In[12]:


def NewScore8_Calc(Prediction_Score,thresh,Quality_Score):
            NewScore8= [] 
            newscore8=0
            b = 0
            for b in range(0,len(Prediction_Score), 1):
                if Prediction_Score[b] > thresh: newscore8= Quality_Score[b] 
                else: newscore8=Prediction_Score[b]
                NewScore8.append(newscore8)
            return NewScore8


# In[13]:


iups_NewScore8=NewScore8_Calc(iups_PredictionScore,0.5,iups_QualityScore)
disopred_NewScore8=NewScore8_Calc(disopred_PredictionScore,0.5,disopred_QualityScore)
spot_NewScore8=NewScore8_Calc(spot_PredictionScore,0.49,spot_QualityScore)


# In[14]:


iups_NewScore8_byProt= np.split(iups_NewScore8,d)
disopred_NewScore8_byProt= np.split(disopred_NewScore8,d)
spot_NewScore8_byProt= np.split(spot_NewScore8,d)


# In[15]:


def Window_Calculation(NewScore10_byProteins):
              Prediction_ScoreWindow= [] 
              newscore10window=0
              b = 0
              for NewScore10_byProtein in NewScore10_byProteins:
                                        for b in range(0, len(NewScore10_byProtein), 1):
                                                      if b == 0 : newscore10window= (NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/6
                                                      elif b == 1 : newscore10window= (NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/7
                                                      elif b == 2 : newscore10window= (NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/8
                                                      elif b == 3 : newscore10window= (NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/9
                                                      elif b == 4 : newscore10window= (NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/10
                                                      elif len(NewScore10_byProtein)-b == 5 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4])/10
                                                      elif len(NewScore10_byProtein)-b == 4 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3])/9
                                                      elif len(NewScore10_byProtein)-b == 3 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2])/8
                                                      elif len(NewScore10_byProtein)-b == 2 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1])/7
                                                      elif len(NewScore10_byProtein)-b == 1 : newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b])/6  
                                                      else: newscore10window= (NewScore10_byProtein[b-5]+NewScore10_byProtein[b-4]+NewScore10_byProtein[b-3]+NewScore10_byProtein[b-2]+NewScore10_byProtein[b-1]+NewScore10_byProtein[b]+NewScore10_byProtein[b+1]+NewScore10_byProtein[b+2]+NewScore10_byProtein[b+3]+NewScore10_byProtein[b+4]+NewScore10_byProtein[b+5])/11
                                                      Prediction_ScoreWindow.append(newscore10window)
              return Prediction_ScoreWindow


# In[16]:


def WeightedWindow_Calculation(NewScore10_byProteins):
              Prediction_ScoreWindowWeighted= [] 
              newscore10windowweighted=0
              b = 0
              for NewScore10_byProtein in NewScore10_byProteins:
                                        for b in range(0, len(NewScore10_byProtein), 1):
                                                      if b == 0 : newscore10windowweighted= (NewScore10_byProtein[b]*0.65+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/6
                                                      elif b == 1 : newscore10windowweighted= (NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.56+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/7
                                                      elif b == 2 : newscore10windowweighted= (NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.48+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/8
                                                      elif b == 3 : newscore10windowweighted= (NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.41+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/9
                                                      elif b == 4 : newscore10windowweighted= (NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.35+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/10
                                                      elif len(NewScore10_byProtein)-b == 5 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.35+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06)/10
                                                      elif len(NewScore10_byProtein)-b == 4 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.41+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07)/9
                                                      elif len(NewScore10_byProtein)-b == 3 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.48+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08)/8
                                                      elif len(NewScore10_byProtein)-b == 2 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.56+NewScore10_byProtein[b+1]*0.09)/7
                                                      elif len(NewScore10_byProtein)-b == 1 : newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.65)/6  
                                                      else: newscore10windowweighted= (NewScore10_byProtein[b-5]*0.05+NewScore10_byProtein[b-4]*0.06+NewScore10_byProtein[b-3]*0.07+NewScore10_byProtein[b-2]*0.08+NewScore10_byProtein[b-1]*0.09+NewScore10_byProtein[b]*0.3+NewScore10_byProtein[b+1]*0.09+NewScore10_byProtein[b+2]*0.08+NewScore10_byProtein[b+3]*0.07+NewScore10_byProtein[b+4]*0.06+NewScore10_byProtein[b+5]*0.05)/11
                                                      Prediction_ScoreWindowWeighted.append(newscore10windowweighted)
              return Prediction_ScoreWindowWeighted


# In[17]:


def ScoreUp_Calculation(Score_byProteins):
          Prediction_Score_Up=[]
          scoreup=0
          b = 0
          for Score_byProtein in Score_byProteins:
                      for b in range(0,len(Score_byProtein),1):
                                     if b == 0: scoreup =  Score_byProtein[0]
                                     else: scoreup =  Score_byProtein[b-1]
                                     Prediction_Score_Up.append(scoreup)
          return Prediction_Score_Up


# In[18]:


def ScoreDown_Calculation(Score_byProteins):
          Prediction_Score_Down=[]
          scoredown=0
          b = 0
          for Score_byProtein in Score_byProteins:
                      for b in range(0,len(Score_byProtein),1):
                                     if len(Score_byProtein)-b == 1 : scoredown =  Score_byProtein[b]
                                     else: scoredown =  Score_byProtein[b+1]
                                     Prediction_Score_Down.append(scoredown)
          return Prediction_Score_Down


# In[19]:


def Score3AVG_Calculation(Score_byProteins):
          Prediction_Score_3AVG =[]
          score3avg = 0
          b = 0
          for Score_byProtein in Score_byProteins:
                      for b in range(0,len(Score_byProtein),1):
                                     if b == 0: score3avg =  (Score_byProtein[b] + Score_byProtein[b+1] + Score_byProtein[b+1])/3
                                     elif len(Score_byProtein)-b == 1 : score3avg =  (Score_byProtein[b] + Score_byProtein[b-1] + Score_byProtein[b-1])/3
                                     else: score3avg =  (Score_byProtein[b-1] + Score_byProtein[b] + Score_byProtein[b+1])/3
                                     Prediction_Score_3AVG.append(score3avg)
          return Prediction_Score_3AVG


# In[20]:


def DisoderDept_Calc(Binary_Prediction):
          Disorder_Index= [] 
          disorderindex=0
          b = 0
          for b in range(0, len(Binary_Prediction), 1):
                                                   if Binary_Prediction[b] == 1 : disorderindex= b
                                                   Disorder_Index.append(disorderindex)
          Disorder_Index=np.unique(Disorder_Index)
          
          ######################################################################################################
          #Getting the margins of disorder regions
          Disorder_Margins= [] 
          disordermargins=0
          b = 0
          h = len(Disorder_Index) -1
          for b in range(0, h , 1):
                                      if Disorder_Index[b+1]-Disorder_Index[b] != 1 : disordermargins= (b+1)
                                      elif np.isin(Disorder_Index[b], d) == True : disordermargins= (b)
                                      Disorder_Margins.append(disordermargins)
          Disorder_Margins=np.unique(Disorder_Margins)
          
          ##############################################################################################################
          #Splitting the disorder positions in to regions
          Disorder_Regions= np.split(Disorder_Index,Disorder_Margins)
          
          ##############################################################################################################
          #Calculating the depth of all disorder residues in their regons
          Disorder_Depth= [] 
          disorderdepth=0
          b = 0
          for Disorder_Region in Disorder_Regions:
                                          for b in range(0, len(Disorder_Region), 1):
                                                                if b < (len(Disorder_Region))/2 : disorderdepth = b 
                                                                elif b > (len(Disorder_Region))/2 : disorderdepth = len(Disorder_Region) -(b+1) 
                                                                Disorder_Depth.append(disorderdepth)
          
          ################################################################################################################
          #Assigning -1 to ordered residues
          Disorder_DepthScore= [] 
          disorderdepthscore=0
          b = 0
          for b in range(0, len(Binary_Prediction), 1):
                                                   if Binary_Prediction[b] == 0 : disorderdepthscore= -1
                                                   else: disorderdepthscore= Binary_Prediction[b]
                                                   Disorder_DepthScore.append(disorderdepthscore)
                                  
          #############################################################################################################
          #Casting lists to arrays 
          Disorder_DepthScore = np.asarray(Disorder_DepthScore)
          Disorder_Index = np.asarray(Disorder_Index)
          Disorder_Depth = np.asarray(Disorder_Depth)
          
          #Inserting depths to the position of disorder residues
          np.put(Disorder_DepthScore, Disorder_Index, Disorder_Depth)
          Disorder_DepthScore=list(Disorder_DepthScore)
          return Disorder_DepthScore


# In[21]:


def OderDept_Calc(Binary_Prediction):
          Order_Index= [] 
          orderindex=0
          b = 0
          for b in range(0, len(Binary_Prediction), 1):
                                                   if Binary_Prediction[b] == 0 : orderindex= b
                                                   Order_Index.append(orderindex)
          Order_Index=np.unique(Order_Index)
          
          ######################################################################################################
          #Getting the margins of disorder regions
          Order_Margins= [] 
          ordermargins=0
          b = 0
          h = len(Order_Index) -1
          for b in range(0,h , 1):
                                      if Order_Index[b+1]-Order_Index[b] != 1 : ordermargins= (b+1)
                                      elif np.isin(Order_Index[b], d) == True : ordermargins= (b)
                                      Order_Margins.append(ordermargins)
          Order_Margins=np.unique(Order_Margins)
          
          #############################################################################################################
          #Splitting the disorder positions in to regions
          Order_Regions= np.split(Order_Index,Order_Margins)
          
          #############################################################################################################
          #Calculating the depth of all disorder residues in their regons
          Order_Depth= [] 
          orderdepth=0
          b = 0
          for Order_Region in Order_Regions:
                                          for b in range(0, len(Order_Region), 1):
                                                                if b < (len(Order_Region))/2 : orderdepth = b 
                                                                elif b > (len(Order_Region))/2 : orderdepth = len(Order_Region) -(b+1)
                                                                Order_Depth.append(orderdepth)
          
          ################################################################################################################
          #Assigning -1 to disordered residues
          Order_DepthScore= [] 
          orderdepthscore=0
          b = 0
          for b in range(0, len(Binary_Prediction), 1):
                                                   if Binary_Prediction[b] == 1 : orderdepthscore= -1
                                                   else: orderdepthscore= Binary_Prediction[b]
                                                   Order_DepthScore.append(orderdepthscore)
                                  
          #############################################################################################################
          #Casting lists to arrays 
          Order_DepthScore = np.asarray(Order_DepthScore)
          Order_Index = np.asarray(Order_Index)
          Order_Depth = np.asarray(Order_Depth)
          
          #Inserting depths to the position of disorder residues
          np.put(Order_DepthScore, Order_Index, Order_Depth)
          Order_DepthScore=list(Order_DepthScore)
          return Order_DepthScore


# In[22]:


def TerminalDistance_Calc(Prediction_Score_byProteins):
     Terminal_Distance= [] 
     terminaldistance=0
     b = 0
     for Prediction_Score_byProtein in Prediction_Score_byProteins:
                                             for b in range(0, (len(Prediction_Score_byProtein)), 1):
                                                           if b <len(Prediction_Score_byProtein)/2 : terminaldistance= b/len(Prediction_Score_byProtein)
                                                           else : terminaldistance = ((len(Prediction_Score_byProtein))-(b+1))/len(Prediction_Score_byProtein) 
                                                           Terminal_Distance.append(terminaldistance)
     return Terminal_Distance


# In[23]:


def TerminalDistance10_Calc(Prediction_Score_byProteins):
      Terminal_Distance10= [] 
      terminaldistance10=0
      b = 0
      for Prediction_Score_byProtein in Prediction_Score_byProteins:
                                        for b in range(0, (len(Prediction_Score_byProtein)), 1):
                                                      if b <10 : terminaldistance10= b
                                                      elif (len(Prediction_Score_byProtein))-(b+1)<10 :terminaldistance10 = ((len(Prediction_Score_byProtein))-(b+1))
                                                      else : terminaldistance10 = 10
                                                      Terminal_Distance10.append(terminaldistance10)
      return Terminal_Distance10


# In[24]:


iups_Prediction_ScoreWindow=Window_Calculation(iups_PredictionScore_byProt)
iups_Prediction_ScoreWindowWeighted=WeightedWindow_Calculation(iups_PredictionScore_byProt)
iups_Prediction_Score_Up=ScoreUp_Calculation(iups_PredictionScore_byProt) 
iups_Prediction_Score_Down= ScoreDown_Calculation(iups_PredictionScore_byProt)
iups_Prediction_Score_3AVG=Score3AVG_Calculation(iups_PredictionScore_byProt)
iups_NewScore8Window=Window_Calculation(iups_NewScore8_byProt) 
iups_NewScore8WindowWeighted=WeightedWindow_Calculation(iups_NewScore8_byProt) 
iups_NewScore8_Up =ScoreUp_Calculation(iups_NewScore8_byProt) 
iups_NewScore8_Down=ScoreDown_Calculation(iups_NewScore8_byProt)
iups_NewScore8_3AVG=Score3AVG_Calculation(iups_NewScore8_byProt)
iups_Disorder_DepthScore=DisoderDept_Calc(iups_Binary)
iups_Order_DepthScore = OderDept_Calc(iups_Binary)
iups_Terminal_Distance= TerminalDistance_Calc(iups_PredictionScore_byProt)
iups_Terminal_Distance10=TerminalDistance10_Calc(iups_PredictionScore_byProt)


# In[25]:


disopred_Prediction_ScoreWindow=Window_Calculation(disopred_PredictionScore_byProt)
disopred_Prediction_ScoreWindowWeighted=WeightedWindow_Calculation(disopred_PredictionScore_byProt)
disopred_Prediction_Score_Up=ScoreUp_Calculation(disopred_PredictionScore_byProt) 
disopred_Prediction_Score_Down= ScoreDown_Calculation(disopred_PredictionScore_byProt)
disopred_Prediction_Score_3AVG=Score3AVG_Calculation(disopred_PredictionScore_byProt)
disopred_NewScore8Window=Window_Calculation(disopred_NewScore8_byProt) 
disopred_NewScore8WindowWeighted=WeightedWindow_Calculation(disopred_NewScore8_byProt) 
disopred_NewScore8_Up =ScoreUp_Calculation(disopred_NewScore8_byProt) 
disopred_NewScore8_Down=ScoreDown_Calculation(disopred_NewScore8_byProt)
disopred_NewScore8_3AVG=Score3AVG_Calculation(disopred_NewScore8_byProt)
disopred_Disorder_DepthScore=DisoderDept_Calc(disopred_Binary)
disopred_Order_DepthScore = OderDept_Calc(disopred_Binary)
disopred_Terminal_Distance= TerminalDistance_Calc(disopred_PredictionScore_byProt)
disopred_Terminal_Distance10=TerminalDistance10_Calc(disopred_PredictionScore_byProt)


# In[26]:


spot_Prediction_ScoreWindow=Window_Calculation(spot_PredictionScore_byProt)
spot_Prediction_ScoreWindowWeighted=WeightedWindow_Calculation(spot_PredictionScore_byProt)
spot_Prediction_Score_Up=ScoreUp_Calculation(spot_PredictionScore_byProt) 
spot_Prediction_Score_Down= ScoreDown_Calculation(spot_PredictionScore_byProt)
spot_Prediction_Score_3AVG=Score3AVG_Calculation(spot_PredictionScore_byProt)
spot_NewScore8Window=Window_Calculation(spot_NewScore8_byProt) 
spot_NewScore8WindowWeighted=WeightedWindow_Calculation(spot_NewScore8_byProt) 
spot_NewScore8_Up =ScoreUp_Calculation(spot_NewScore8_byProt) 
spot_NewScore8_Down=ScoreDown_Calculation(spot_NewScore8_byProt)
spot_NewScore8_3AVG=Score3AVG_Calculation(spot_NewScore8_byProt)
spot_Disorder_DepthScore=DisoderDept_Calc(spot_Binary)
spot_Order_DepthScore = OderDept_Calc(spot_Binary)
spot_Terminal_Distance= TerminalDistance_Calc(spot_PredictionScore_byProt)
spot_Terminal_Distance10=TerminalDistance10_Calc(spot_PredictionScore_byProt)


# In[27]:


Featuer_Frame=pd.DataFrame({
'iups_PredictionScore':iups_PredictionScore,
'iups_Prediction_ScoreWindow':iups_Prediction_ScoreWindow,
'iups_Prediction_ScoreWindowWeighted':iups_Prediction_ScoreWindowWeighted,
'iups_Prediction_Score_Up':iups_Prediction_Score_Up,
'iups_Prediction_Score_Down':iups_Prediction_Score_Down,
'iups_Prediction_Score_3AVG':iups_Prediction_Score_3AVG,
'iups_NewScore8':iups_NewScore8,
'iups_NewScore8Window':iups_NewScore8Window,
'iups_NewScore8WindowWeighted':iups_NewScore8WindowWeighted,
'iups_NewScore8_Up':iups_NewScore8_Up, 
'iups_NewScore8_Down':iups_NewScore8_Down,
'iups_NewScore8_3AVG':iups_NewScore8_3AVG,
'iups_Disorder_DepthScore':iups_Disorder_DepthScore,
'iups_Order_DepthScore':iups_Order_DepthScore, 
'iups_Terminal_Distance':iups_Terminal_Distance,
'iups_Terminal_Distance10':iups_Terminal_Distance10,
'spot_PredictionScore':spot_PredictionScore,
'spot_Prediction_ScoreWindow':spot_Prediction_ScoreWindow,
'spot_Prediction_ScoreWindowWeighted':spot_Prediction_ScoreWindowWeighted,
'spot_Prediction_Score_Up':spot_Prediction_Score_Up,
'spot_Prediction_Score_Down':spot_Prediction_Score_Down,
'spot_Prediction_Score_3AVG':spot_Prediction_Score_3AVG,
'spot_NewScore8':spot_NewScore8,
'spot_NewScore8Window':spot_NewScore8Window,
'spot_NewScore8WindowWeighted':spot_NewScore8WindowWeighted,
'spot_NewScore8_Up':spot_NewScore8_Up, 
'spot_NewScore8_Down':spot_NewScore8_Down,
'spot_NewScore8_3AVG':spot_NewScore8_3AVG,
'spot_Disorder_DepthScore':spot_Disorder_DepthScore,
'spot_Order_DepthScore':spot_Order_DepthScore, 
'spot_Terminal_Distance':spot_Terminal_Distance,
'spot_Terminal_Distance10':spot_Terminal_Distance10,
'disopred_PredictionScore':disopred_PredictionScore,
'disopred_Prediction_ScoreWindow':disopred_Prediction_ScoreWindow,
'disopred_Prediction_ScoreWindowWeighted':disopred_Prediction_ScoreWindowWeighted,
'disopred_Prediction_Score_Up':disopred_Prediction_Score_Up,
'disopred_Prediction_Score_Down':disopred_Prediction_Score_Down,
'disopred_Prediction_Score_3AVG':disopred_Prediction_Score_3AVG,
'disopred_NewScore8':disopred_NewScore8,
'disopred_NewScore8Window':disopred_NewScore8Window,
'disopred_NewScore8WindowWeighted':disopred_NewScore8WindowWeighted,
'disopred_NewScore8_Up':disopred_NewScore8_Up, 
'disopred_NewScore8_Down':disopred_NewScore8_Down,
'disopred_NewScore8_3AVG':disopred_NewScore8_3AVG,
'disopred_Disorder_DepthScore':disopred_Disorder_DepthScore,
'disopred_Order_DepthScore':disopred_Order_DepthScore, 
'disopred_Terminal_Distance':disopred_Terminal_Distance,
'disopred_Terminal_Distance10':disopred_Terminal_Distance10,
'iups_mcc_Predicted':iups_mcc_Predicted,
'spot_mcc_Predicted':spot_mcc_Predicted,
'disopred_mcc_Predicted':disopred_mcc_Predicted})
Featuer_Frame.head()


# In[28]:


X_test = np.array(Featuer_Frame.iloc[:,0:51])
Xtest_3d = X_test.reshape(len(disopred_mcc_Predicted), 51, 1)


# In[29]:


model = tf.keras.models.load_model('LSTM_MCC.h5')


# In[30]:

with tf.device('/cpu:0'):
       quarter_PredictionScore = model.predict(Xtest_3d, batch_size=100)
quarter_PredictionScore=np.concatenate(quarter_PredictionScore)


# In[31]:


quarter_Binary=Score_Binary(quarter_PredictionScore,0.4600)
quarter_Binary_byProt=np.split(quarter_Binary,d)
quarter_PredictionScore_byProt=np.split(quarter_PredictionScore,d)


# In[32]:


dis465_MCC_byProt=np.split(dis465_mcc_Predicted,d)
disHL_MCC_byProt=np.split(disHL_mcc_Predicted,d)
espD_MCC_byProt=np.split(espD_mcc_Predicted,d)
espN_MCC_byProt=np.split(espN_mcc_Predicted,d)
espX_MCC_byProt=np.split(espX_mcc_Predicted,d)
glo_MCC_byProt=np.split(glo_mcc_Predicted,d)
iupl_MCC_byProt=np.split(iupl_mcc_Predicted,d)
iups_MCC_byProt=np.split(iups_mcc_Predicted,d)
vsl_MCC_byProt=np.split(vsl_mcc_Predicted,d)
spot_MCC_byProt=np.split(spot_mcc_Predicted,d)
disopred_MCC_byProt=np.split(disopred_mcc_Predicted,d)


# In[33]:


Lines=[]
b=0
for b in range(0,len(Uniprot_ID),1):
                Lines.append(Uniprot_ID[b])
                Lines.append(AA_Sequence[b])
                x = str(list(dis465_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(dis465_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(dis465_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(disHL_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(disHL_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(disHL_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(espD_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(espD_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(espD_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(espN_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(espN_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(espN_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(espX_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(espX_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(espX_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(glo_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(glo_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(glo_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(iupl_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(iupl_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(iupl_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(iups_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(iups_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(iups_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(vsl_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(vsl_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(vsl_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                
                x = str(list(spot_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(spot_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(spot_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(disopred_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(disopred_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(disopred_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)
                
                x = str(list(quarter_Binary_byProt[b]))
                x = x[1:-1]
                x = x.replace(",", "")
                x = x.replace(" ", "")
                Lines.append(x)
                x = str(list(np.around(quarter_PredictionScore_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)


# In[34]:


with open('Qplus_Output.txt', 'w') as f:
    for item in Lines:
        f.write("%s\n" % item)

