# -*- coding=utf-8 -*-
import os
import sys
import pandas as pd
import numpy as np 

patientInfo=pd.read_csv('clinical_info',sep='\t')
patientInfo=patientInfo.loc[patientInfo['Biopsy_Time']=='pre-treatment']

def countnsmut():
    EmutNumber= dict()
    for id in patientInfo['Run']:
        if(os.path.isfile("D:/ICB_data/somatic_nonsynonymous_mut/"+id+".read_count")):
        #Checking to see if the file exis  
            dtmp=pd.read_csv("D:/ICB_data/somatic_nonsynonymous_mut/"+id+".read_count",sep='\t')
            dtmp=dtmp.drop_duplicates()  #去重复
            EmutNumber[id]=(len(dtmp[dtmp["alt_count"]>4]))
            dtmp.drop(dtmp.index,inplace=True)
    return (EmutNumber)



MutNumber=countnsmut()
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
patientInfo.to_excel("./join_nun_info.xlsx",index=False)


#separate cance type 、anti target and response type

patientInfo[(patientInfo.Cancer=="lung cancer") & (patientInfo.Anti_target=="anti-PD1") ].to_csv("LungCancer.csv",index=False)
patientInfo[(patientInfo.Cancer=="melanoma") & (patientInfo.Anti_target=="anti-CTLA4") ].to_csv("Melanoma_CTLA4.csv",index=False)
patientInfo[(patientInfo.Cancer=="melanoma") & (patientInfo.Anti_target=="anti-PD1") ].to_csv("Melanoma_PD1.csv",index=False)


LungCancer=pd.read_csv("LungCancer.csv",sep=',')
Melanoma_CTLA4=pd.read_csv("Melanoma_CTLA4.csv",sep=',')
Melanoma_PD1=pd.read_csv("Melanoma_PD1.csv",sep=',')



#去掉离群值：肺癌样本中去掉突变数目大于400的样本，melanoma-pd1去掉突变数目大于1000的样本，melanoma CTLA4去掉大于2000的样本。

LungCancer[["All_mutations"]]=LungCancer[["All_mutations"]].astype('int')  #强制类型转换
LungCancer=LungCancer[LungCancer.All_mutations<=400]

Melanoma_PD1[["All_mutations"]]=LungCancer[["All_mutations"]].astype('int')
Melanoma_PD1=Melanoma_PD1[Melanoma_PD1.All_mutations<=1000]

Melanoma_CTLA4[["All_mutations"]]=LungCancer[["All_mutations"]].astype('int')
Melanoma_CTLA4=Melanoma_CTLA4[Melanoma_CTLA4.All_mutations<=2000]



"""

box1=LungCancer[['Run']][LungCancer.Response_type=='PR']
box2=LungCancer[['Run']][LungCancer.Response_type=='PD']
box1['Response_type']='CR/PR'   #Response type for box1 is PR
box2['Response_type']='PD'
Btmp=box1.append(box2)
Btmp=Btmp.reset_index(drop=True)#重新设置列索引



response_sign="CR/PR"
filepath="Melanoma_PD1_mutation_gene_matrix"
box1=Melanoma_PD1[['Run']][(Melanoma_PD1.Response_type=='CR') | (Melanoma_PD1.Response_type=='PR') ]
box2=Melanoma_PD1[['Run']][Melanoma_PD1.Response_type=='PD' ]
box1['Response_type']='CR/PR'   
box2['Response_type']='PD'
Btmp=box1.append(box2)

box1=Melanoma_CTLA4[['Run']][Melanoma_CTLA4.Response_type=='long-term-benefit']
box2=Melanoma_CTLA4[['Run']][Melanoma_CTLA4.Response_type=='minimal-or-no-benefit']
box1['Response_type']='long-term-benefit'
box2['Response_type']='minimal-or-no-benefit'
Btmp=box1.append(box2)
Btmp=Btmp.reset_index(drop=True)#重新设置列索引

"""

sample_mutation_gene_dict={}
all_mutation_genes=[]
mutation_gene_list=[]

response_sign="long-term-benefit"
filepath="Melanoma_CTLA4_mutation_gene_matrix"
box1=Melanoma_CTLA4[['Run']][Melanoma_CTLA4.Response_type=='long-term-benefit']
box2=Melanoma_CTLA4[['Run']][Melanoma_CTLA4.Response_type=='minimal-or-no-benefit']
box1['Response_type']='long-term-benefit'
box2['Response_type']='minimal-or-no-benefit'
Btmp=box1.append(box2)
Btmp=Btmp.reset_index(drop=True)#重新设置列索引



Btmp=Btmp.reset_index(drop=True)#重新设置列索引



for id in Btmp["Run"]:
    mutation_gene=pd.read_csv("D:/ICB_data/somatic_nonsynonymous_mut/"+id+"_convert_vcf.txt",sep='\t')
    mutation_gene_list=mutation_gene["ensembl_gene_id"].drop_duplicates().values.tolist()
    sample_mutation_gene_dict[id]=mutation_gene_list
    all_mutation_genes=list(set(all_mutation_genes).union(set(mutation_gene_list)))
    
    mutation_gene_list=[]
    mutation_gene.drop(mutation_gene.index,inplace=True)

all_mutation_genes.append("response")
mutation_gene_matrix=pd.DataFrame(columns=all_mutation_genes)

for s_id,i in zip(Btmp["Run"],range(len(Btmp["Run"]))):
    mutation_gene_list=sample_mutation_gene_dict[s_id]
    if Btmp.loc[i,"Response_type"]==response_sign:
        s_response=1
    else:
        s_response=0
    list=np.ones((1,len(mutation_gene_list)),dtype=int)
    sample_mutation_gene_matrix=pd.DataFrame(list,columns=mutation_gene_list)
    sample_mutation_gene_matrix["response"]=s_response
    mutation_gene_matrix=pd.concat([mutation_gene_matrix,sample_mutation_gene_matrix],ignore_index=True)
    sample_mutation_gene_matrix.drop(sample_mutation_gene_matrix.index,inplace=True)
    mutation_gene_list=[]

mutation_gene_matrix.fillna(0,inplace=True)
print(mutation_gene_matrix)
mutation_gene_matrix.to_csv(filepath+".csv",index=False)
print(filepath+" write to file!")




