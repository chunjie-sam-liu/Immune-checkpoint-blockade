# -*- coding=utf-8 -*-
import os
import subprocess
import pandas as pd


patientInfo=pd.read_csv('clinical_info',sep='\t')
patientInfo=patientInfo.loc[patientInfo['Biopsy_Time']=='pre-treatment']

def countmut():
    #Count the number of somatic mutations
    MutNumber = dict()
    EmutNumber= dict()
    for id in patientInfo['Run']:
        if(os.path.isfile("./somatic_mut/"+id+".infoForPyclone")):
        #Checking to see if the file exis
            dtmp=pd.read_csv("./somatic_mut/"+id+".infoForPyclone",sep='\t',header=None)
            MutNumber[id]=(len(dtmp))
            EmutNumber[id]=(len(dtmp[dtmp[6]>4]))
            dtmp.drop(dtmp.index,inplace=True)
    return (MutNumber,EmutNumber)

def countnsmut():
    MutNumber = dict()
    EmutNumber= dict()
    for id in patientInfo['Run']:
        if(os.path.isfile("./somatic_nonsynonymous_mut/"+id+".read_count")):
        #Checking to see if the file exis  
            dtmp=pd.read_csv("./somatic_nonsynonymous_mut/"+id+".read_count",sep='\t')
            dtmp=dtmp.drop_duplicates()  #去重复
            MutNumber[id]=(len(dtmp))
            EmutNumber[id]=(len(dtmp[dtmp["alt_count"]>4]))
            dtmp.drop(dtmp.index,inplace=True)
    return (MutNumber,EmutNumber)


Mut=dict()
Emut=dict()
muttype=2         
if (muttype==1):
    (Mut,Emut)=countmut()
    print("体细胞突变")
else:
    (Mut,Emut)=countnsmut()
    print("非同义突变")

MutNumber=dict()

Type=2
if (Type==1):
    MutNumber=Mut
    print("不同反应群体突变数的差异 ")
else:
    MutNumber=Emut
    print("表达的突变（第七列的值大于4）在不同反应群体之间的差异")
    

#Add the number of mutations to the patientinfo sheet
stmp1 = pd.Series(MutNumber)
dict_tmp1 = {'Run':stmp1.index,'All_mutations':stmp1.values}
df_tmp1 = pd.DataFrame(dict_tmp1)
patientInfo=pd.merge(patientInfo,df_tmp1,left_on='Run',right_on='Run')

#Determine the type of response add this to clinical_info
ResponseType=dict()
for i in range(0, len(patientInfo)-1):
    if(pd.isnull(patientInfo.iloc[i]['Response']) and pd.isnull(patientInfo.iloc[i]['Second_Response_standard'])):
        patientInfo.drop(i)
    elif(pd.notnull(patientInfo.iloc[i]['Response'])):
        ResponseType[patientInfo.iloc[i]['Run']]=patientInfo.iloc[i]['Response']
    elif(pd.isnull(patientInfo.iloc[i]['Response'])  and pd.notnull(patientInfo.iloc[i]['Second_Response_standard'])):
        ResponseType[patientInfo.iloc[i]['Run']]=patientInfo.iloc[i]['Second_Response_standard']
stmp2 = pd.Series(ResponseType)
dict_tmp2 = {'Run':stmp2.index,'Response_type':stmp2.values}
df_tmp2 = pd.DataFrame(dict_tmp2)
patientInfo=pd.merge(patientInfo,df_tmp2,left_on='Run',right_on='Run')
patientInfo.to_excel("./tmp/join_nun_info.xlsx")


#separate cance type 、anti target and response type
Filter = dict(list(patientInfo.groupby(['Cancer','Anti_target'])))
LungCancer=Filter['lung cancer','anti-PD1'] 
Melanoma_CTLA4=Filter['melanoma','anti-CTLA4']
Melanoma_PD1=Filter['melanoma','anti-PD1']


def lungcancer_data_prepare_drawboxplot():
    box1=LungCancer[['All_mutations']][LungCancer.Response_type=='PR']
    box2=LungCancer[['All_mutations']][LungCancer.Response_type=='PD']
    box3=LungCancer[['All_mutations']][LungCancer.Response_type=='SD']
    #Set CR and PR as CR/PR
    box1['Response_type']='CR/PR'   #Response type for box1 is PR
    box2['Response_type']='PD'
    box3['Response_type']='SD'
    Btmp=box1.append(box2)
    Btmp=Btmp.append(box3)
    Btmp.to_csv('./tmp/lung_cancer.csv')


def Melanoma_CTLA4_data_prepare_drawboxplot():
    box1=Melanoma_CTLA4[['All_mutations']][Melanoma_CTLA4.Response_type=='long-term-benefit']
    box2=Melanoma_CTLA4[['All_mutations']][Melanoma_CTLA4.Response_type=='minimal-or-no-benefit']
    box1['Response_type']='long-term-benefit'
    box2['Response_type']='minimal-or-no-benefit'
    box1=box1.append(box2)
    box1.to_csv('./tmp/Melanoma_CTLA4.csv')



def Melanoma_PD1_data_prepare_drawboxplot():
    box1=Melanoma_PD1[['All_mutations']][(Melanoma_PD1.Response_type=='CR') |(Melanoma_PD1.Response_type=='PR') ]
    box2=Melanoma_PD1[['All_mutations']][Melanoma_PD1.Response_type=='PD' ]
    box1['Response_type']='CR/PR'
    box2['Response_type']='PD'
    box1=box1.append(box2)
    box1.to_csv('./tmp/Melanoma_PD1.csv')
 

def main():
    lungcancer_data_prepare_drawboxplot()
    Melanoma_CTLA4_data_prepare_drawboxplot()
    Melanoma_PD1_data_prepare_drawboxplot()
    print("Done!")
    #Rscript [2]Draw_boxplot.R 1
 
if __name__ == '__main__':
    main()
