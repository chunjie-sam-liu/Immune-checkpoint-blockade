import os
import pandas as pd

patientInfo=pd.read_csv('clinical_info',sep='\t')
patientInfo=patientInfo.loc[patientInfo['Biopsy_Time']=='pre-treatment']

#Count the number of somatic mutations
MutNumber = dict()
for id in patientInfo['Run']:
    if(os.path.isfile("./somatic_mut/"+id+".infoForPyclone")):
    #Checking to see if the file exis
        count=0
        for index, line in enumerate(open("./somatic_mut/"+id+".infoForPyclone",'r')):
        #Calculate number of rows
            count += 1
        MutNumber[id] = count

#Add the number of mutations to the patientinfo sheet 
stmp1 = pd.Series(MutNumber)
dict_tmp1 = {'Run':stmp1.index,'MutNumber':stmp1.values}
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
patientInfo.to_excel("join_nun_info.xlsx")


#根据不同的cance type 、anti target和 response type 区分样本信息
Filter = dict(list(patientInfo.groupby(['Cancer','Anti_target','Response_type'])))


LungCancer_CR_PR=Filter['lung cancer','anti-PD1','PR'] 
LungCancer_PD=Filter['lung cancer','anti-PD1','PD']
LungCancer_SD=Filter['lung cancer','anti-PD1','SD']
Melanoma_CTLA4_Ltb=Filter['melanoma','anti-CTLA4','long-term-benefit']
Melanoma_CTLA4_Monb=Filter['melanoma','anti-CTLA4','minimal-or-no-benefit']
Melanoma_PD1_PR=Filter['melanoma','anti-PD1','PR']
Melanoma_PD1_CR=Filter['melanoma','anti-PD1','CR']
Melanoma_PD1_PD=Filter['melanoma','anti-PD1','PD']
Melanoma_PD1_CR_PR=pd.concat([Melanoma_PD1_PR,Melanoma_PD1_CR])
