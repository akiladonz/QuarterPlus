#!/home/pyenv/versions/py3.7/bin/python
# coding: utf-8

# In[1]:
import numpy as np
import pandas as pd
import bisect
import sys
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os


input_txt_path =  os.path.abspath(sys.argv[1])
output_figure_path = os.path.abspath(sys.argv[2])
output_text_path = os.path.abspath(sys.argv[3])


Lines= []
f = open(input_txt_path,'r')
for line in f:
    data_line = line.rstrip().split('\t')
    Lines.append(data_line)		
	
	
	
Uniprot_ID=[]
AA_Sequence=[]
iups_PredictionScore_byProt=[]
iups_QualityScore_byProt=[]
spot_PredictionScore_byProt=[]
spot_QualityScore_byProt=[]
disopred_PredictionScore_byProt=[]
disopred_QualityScore_byProt=[]
quarter_PredictionScore_byProt=[]
quarter_QualityScore_byProt=[]
J = len(Lines)-9
b =0
for b in range(0,J,10):
    Uniprot_ID.append(Lines[b][0])
    AA_Sequence.append(Lines[b+1][0])
    iups_PredictionScore_byProt.append(list(map(float,Lines[b+2][0].rstrip().split(','))))
    iups_QualityScore_byProt.append(list(map(float,Lines[b+3][0].rstrip().split(','))))
    spot_PredictionScore_byProt.append(list(map(float,Lines[b+4][0].rstrip().split(','))))
    spot_QualityScore_byProt.append(list(map(float,Lines[b+5][0].rstrip().split(','))))
    disopred_PredictionScore_byProt.append(list(map(float,Lines[b+6][0].rstrip().split(','))))
    disopred_QualityScore_byProt.append(list(map(float,Lines[b+7][0].rstrip().split(','))))
    quarter_PredictionScore_byProt.append(list(map(float,Lines[b+8][0].rstrip().split(','))))
    quarter_QualityScore_byProt.append(list(map(float,Lines[b+9][0].rstrip().split(','))))


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



iups_PredictionScore=np.concatenate(iups_PredictionScore_byProt)
iups_QualityScore=np.concatenate(iups_QualityScore_byProt)
spot_PredictionScore=np.concatenate(spot_PredictionScore_byProt)
spot_QualityScore=np.concatenate(spot_QualityScore_byProt)
disopred_PredictionScore=np.concatenate(disopred_PredictionScore_byProt)
disopred_QualityScore=np.concatenate(disopred_QualityScore_byProt)
quarter_PredictionScore=np.concatenate(quarter_PredictionScore_byProt)
quarter_QualityScore=np.concatenate(quarter_QualityScore_byProt)



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

iups_mcc_Predicted=[]
spot_mcc_Predicted=[]
disopred_mcc_Predicted=[]
b=0
for b in range(0,len(disopred_QualityScore),1):

                 try:
                   iups_mcc=iups_MCCWindow_tr[bisect.bisect_left(iups_QAWindow_tr, iups_QualityScore[b])]
                 except:
                   iups_mcc=iups_MCCWindow_tr[(bisect.bisect_left(iups_QAWindow_tr, iups_QualityScore[b]))-1]

                 try:
                   spot_mcc=spot_MCCWindow_tr[bisect.bisect_left(spot_QAWindow_tr, spot_QualityScore[b])]
                 except:
                    spot_mcc=spot_MCCWindow_tr[(bisect.bisect_left(spot_QAWindow_tr, spot_QualityScore[b]))-1]
                 try:
                   disopred_mcc=disopred_MCCWindow_tr[bisect.bisect_left(disopred_QAWindow_tr, disopred_QualityScore[b])]
                 except:
                    disopred_mcc=disopred_MCCWindow_tr[(bisect.bisect_left(disopred_QAWindow_tr, disopred_QualityScore[b]))-1]

                 iups_mcc_Predicted.append(iups_mcc)
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
quarter_Binary=Score_Binary(quarter_PredictionScore,0.46)



# In[11]:


disopred_Binary_byProt= np.split(disopred_Binary,d)
spot_Binary_byProt= np.split(spot_Binary,d)
iups_Binary_byProt= np.split(iups_Binary,d)
quarter_Binary_byProt= np.split(quarter_Binary,d)

	
	
		




file = (r"LSTM_Lookup.csv")
LSTM_Lookup= pd.read_csv(file)	
LSTM_QAWindow_tr=LSTM_Lookup.iloc[:,(1)].values
LSTM_MCCWindow_tr=LSTM_Lookup.iloc[:,(2)].values

quarter_mcc_Predicted=[]
b=0
for b in range(0,len(quarter_QualityScore),1):
                 try:
                   mcc=LSTM_MCCWindow_tr[bisect.bisect_left(LSTM_QAWindow_tr, quarter_QualityScore[b])]
                 except:
                   mcc=LSTM_MCCWindow_tr[(bisect.bisect_left(LSTM_QAWindow_tr, quarter_QualityScore[b]))-1]
                 quarter_mcc_Predicted.append(mcc)



quarter_MCC_byProt=np.split(quarter_mcc_Predicted,d)
spot_MCC_byProt=np.split(spot_mcc_Predicted,d)
disopred_MCC_byProt=np.split(disopred_mcc_Predicted,d)
iups_MCC_byProt=np.split(iups_mcc_Predicted,d)




Position =np.arange(0,len(quarter_PredictionScore_byProt[0]),1)
size=5
fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(211)
ax1.plot(Position, quarter_PredictionScore_byProt[0],'#be03fd',label="Quarter+",linewidth=2)
ax1.plot(Position, disopred_PredictionScore_byProt[0],'coral',label="Disopred3",linewidth=2)
ax1.plot(Position, spot_PredictionScore_byProt[0],'#886806',label="SPOT-DISORDER",linewidth=2)
ax1.plot(Position, iups_PredictionScore_byProt[0],'#0652ff',label="IUPred-short",linewidth=2)
ax1.hlines(y=0.5, xmin=0, xmax=len(Position)-1, linewidth=1.5, color='coral',linestyle ='--')
ax1.hlines(y=0.49, xmin=0, xmax=len(Position)-1, linewidth=1.5, color='#886806',linestyle ='--')
ax1.hlines(y=0.5, xmin=0, xmax=len(Position)-1, linewidth=1.5, color='#0652ff',linestyle ='--')
ax1.hlines(y=0.46, xmin=0, xmax=len(Position)-1, linewidth=1.5, color='#be03fd',linestyle ='--')

#ax.legend(bbox_to_anchor=(1.05, 1), loc=4, borderaxespad=0.)
ax1.grid()
#ax1.set_ylabel('  Quater+,Disopred3,SPOT-DISORDER\n and IUPred-short Predicted MCC', rotation=90, ha="right")
fig.text(0.06,0.75, "Quater+,Disopred3,SPOT-DISORDER\n and IUPred-short Propensity Score", ha="center", va="center", rotation=90, fontsize=10)
ax1.set_ylim([0,1.0])
ax1.set_xlim([0,len(Position)-1])
ax1.tick_params(top=False, bottom=False, left=False, right=True, labelleft=True, labelbottom=False)
ax1.yaxis.set_label_coords(-0.05,1.1)
##################################################################################################################
ax2 = fig.add_subplot(212)
A = ['IUPred-short ExpectedMCC',
     'IUPred-short Prediction',
     'DISOPRED3 ExpectedMCC',
     'DISOPRED3 Prediction',
     'SPOT-DISORDER ExpectedMCC',
     'SPOT-DISORDER Prediction',
     'Quarter+ ExpectedMCC',
     'Quarter+ Prediction']
ordercolor='#a8a495'
disordercolor='black'
cmap = matplotlib.cm.get_cmap('RdYlGn')
#################################################################################################
x= np.concatenate(np.argwhere(quarter_Binary_byProt[0]<1))
y = np.full(len(x), 6.0)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=ordercolor)
#********************************************************************************************
x= np.concatenate(np.argwhere(quarter_Binary_byProt[0]>0))
y = np.full(len(x), 6.0)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=disordercolor)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x =np.arange(0,len(Position),1)
y = np.full(len(Position), 5.5)

source_array=quarter_MCC_byProt[0]
b=0
for b in range(0,len(source_array),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=cmap(source_array[b]))
ax2.CLim = [-1,1];

##################################################################################################
x= np.concatenate(np.argwhere(spot_Binary_byProt[0]<1))
y = np.full(len(x), 4.5)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=ordercolor)
#********************************************************************************************
x= np.concatenate(np.argwhere(spot_Binary_byProt[0]>0))
y = np.full(len(x), 4.5)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=disordercolor)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x =np.arange(0,len(Position),1)
y = np.full(len(Position), 4.0)
source_array=spot_MCC_byProt[0]
b=0
for b in range(0,len(source_array),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=cmap(source_array[b]))
ax2.CLim = [-1,1];

#############################################################################################
x= np.concatenate(np.argwhere(disopred_Binary_byProt[0]<1))
y = np.full(len(x), 3.0)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=ordercolor)
#********************************************************************************************
x= np.concatenate(np.argwhere(disopred_Binary_byProt[0]>0))
y = np.full(len(x), 3.0)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=disordercolor)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x =np.arange(0,len(Position),1)
y = np.full(len(Position), 2.5)
source_array=disopred_MCC_byProt[0]
b=0
for b in range(0,len(source_array),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=cmap(source_array[b]))
ax2.CLim = [-1,1];
#############################################################################################
x= np.concatenate(np.argwhere(iups_Binary_byProt[0]<1))
y = np.full(len(x), 1.5)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=ordercolor)
#********************************************************************************************
x= np.concatenate(np.argwhere(iups_Binary_byProt[0]>0))
y = np.full(len(x), 1.5)
b=0
for b in range(0,len(x),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=disordercolor)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x =np.arange(0,len(Position),1)
y = np.full(len(Position), 1.0)
source_array=iups_MCC_byProt[0]
b=0
for b in range(0,len(source_array),1):
               ax2.plot(x[b], y[b], marker='s', markersize=size, color=cmap(source_array[b]))
ax2.CLim = [-1,1];
#############################################################################################
ax2.set_xlim([0,len(Position)-1])
ax2.yaxis.grid()
ax2.xaxis.grid()
#ax2.set_yticks(np.arange(1,8, 1.0))
ticks=[1.0, 1.5, 2.5, 3.0, 4.0, 4.5, 5.5, 6.0]
ax2.set_yticks(ticks)
ax2.set_yticklabels(A)
ax2.get_yticklabels()[7].set_color("#be03fd")
ax2.get_yticklabels()[6].set_color("#be03fd")
ax2.get_yticklabels()[5].set_color("#886806")
ax2.get_yticklabels()[4].set_color("#886806")
ax2.get_yticklabels()[3].set_color("coral")
ax2.get_yticklabels()[2].set_color("coral")
ax2.get_yticklabels()[1].set_color("#0652ff")
ax2.get_yticklabels()[0].set_color("#0652ff")
ax2.tick_params(top=False, bottom=True, left=False, right=False, labelleft=True, labelbottom= True)
ax2.spines["top"].set_visible(False)
ax1.legend(bbox_to_anchor=(0.01, 1.05), loc=3,ncol=1, borderaxespad=0., fontsize='medium')
legend_elements = [Line2D([0], [0], marker='s', color=disordercolor, label='Disordered Regions',markerfacecolor=disordercolor, markersize=5),
                   Line2D([0], [0], marker='s', color=ordercolor, label='Ordered Regions',markerfacecolor=ordercolor, markersize=5),
                  Line2D([0], [0], marker='_', color='#be03fd', label='Quarter+',markerfacecolor='#be03fd', markersize=2),
                  Line2D([0], [0], marker='_', color='coral', label='Disopred3',markerfacecolor='coral', markersize=2),
                  Line2D([0], [0], marker='_', color='#886806', label='SPOT-DISORDER',markerfacecolor='#886806', markersize=2),
                  Line2D([0], [0], marker='_', color='#0652ff', label='IUPred-short',markerfacecolor='#0652ff', markersize=2)]
ax1.legend(handles=legend_elements, bbox_to_anchor=(0.01, 1.05), loc=3,ncol=1, borderaxespad=0., fontsize='medium')                    
norm = mpl.colors.Normalize(vmin=-1, vmax=1)
cax = plt.axes([0.91, 0.2, 0.020, 0.3])
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_label('Expected MCC')
ax2.set_xlabel('Sequence Position')
fig.subplots_adjust(wspace=0, hspace=0.01)
plt.savefig(output_figure_path, bbox_inches='tight',quality =95,orientation ='landscape',dpi =800)









Lines=[]
b=0
for b in range(0,len(Uniprot_ID),1):
                Lines.append(Uniprot_ID[b])
                Lines.append(AA_Sequence[b])

                
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
                x = str(list(np.around(quarter_MCC_byProt[b], decimals=3)))
                x = x[1:-1]
                x = x.replace(" ", "")
                Lines.append(x)


# In[34]:


with open(output_text_path, 'w') as f:
    for item in Lines:
        f.write("%s\n" % item)
