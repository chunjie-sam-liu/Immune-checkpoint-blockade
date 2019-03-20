import os
import re
import glob
import copy
import numpy as np
import pandas as pd 
from collections import Counter
from ShowProcess import ShowProcess

patientInfo=pd.read_csv('clinical_info',sep='\t',skipinitialspace=True )
patientInfo=patientInfo.loc[patientInfo['Biopsy_Time']=='pre-treatment']
sample_list=pd.read_csv('./neoantigen_filter_SB_result/sample_list')
sample_list = sample_list.drop_duplicates()  #去重复

patientInfo=pd.merge(patientInfo,sample_list,left_on='Run',right_on='sample') 
patientInfo.drop(columns=['sample'],inplace=True)

def all_np(arr): #统计list 元素重复的个数。返回一个字典
    arr = np.array(arr)
    key = np.unique(arr)
    result = {}
    for k in key:
        mask = (arr == k)
        arr_new = arr[mask]
        v = arr_new.size
        result[k] = v
    return (result)

hla_list=[]
neoantigen_list=[]
protein_source_list=[]
hla_dict={}
neoantigen_dict={}
protein_source_dict={}
filelist=(glob.glob(r'./neoantigen_filter_SB_result/*.filter'))

max_steps=len(filelist)
process_bar = ShowProcess(max_steps, 'OK')
print("统计基因、新抗原、新抗原蛋白质信息：")
#获得相应的sample ID 与之对应的基因等信息，并计数。
for f in filelist:
    process_bar.show_process()  #打印进度条
    sample_id=(f.split('_')[3].split('\\')[1])  #在Linux中注意斜杠的方向
    if (sample_id in patientInfo['Run'].tolist() and os.path.getsize(f)!=0):
        result_tmp=open(f)
        for line in result_tmp.readlines():
            line = re.sub(' +', '#', line)  #解决分隔符不统一的问题，将出现空格的地方替换为 “#”
            line=line.split("#")
            hla_list.append(line[19])  #添加基因信息
            neoantigen_list.append(line[20]) #添加新抗原信息
            protein_source_list.append(line[28])#添加新的蛋白质的信息

        if  sample_id in hla_dict.keys():
            temp1 = all_np(hla_list)
            temp2 = hla_dict[sample_id]    #这是一个字典  充分利用python的浅拷贝原理改变字典的值
            for k in np.unique(hla_list):
                if (k not in temp2.keys()):
                    temp2[k] = temp1[k]
                temp2[k] += temp1[k]
            hla_list.clear()
        else:
            hla_dict[sample_id]=all_np(hla_list)
            hla_list.clear()

        if  sample_id in neoantigen_dict.keys():
            temp1 = all_np(neoantigen_list)
            temp2 = neoantigen_dict[sample_id]   
            for k in np.unique(neoantigen_list):
                if (k not in temp2.keys()):
                    temp2[k] = temp1[k]
                temp2[k] += temp1[k]
            neoantigen_list.clear()
        else:
            neoantigen_dict[sample_id]=all_np(neoantigen_list)
            neoantigen_list.clear()

        if  sample_id in protein_source_dict.keys():
            temp1 = all_np(protein_source_list)
            temp2 = protein_source_dict[sample_id]   
            for k in np.unique(protein_source_list):
                if (k not in temp2.keys()):
                    temp2[k] = temp1[k]
                temp2[k] += temp1[k]  
            neoantigen_list.clear()
        else:
            protein_source_dict[sample_id]=all_np(protein_source_list)
            protein_source_list.clear()
    
#Determine the type of response add this to clinical_info
ResponseType=dict()
for i in range(0, len(patientInfo)-1):
    if(pd.notnull(patientInfo.iloc[i]['Response'])):
        ResponseType[patientInfo.iloc[i]['Run']]=patientInfo.iloc[i]['Response']
    elif(pd.isnull(patientInfo.iloc[i]['Response'])  and pd.notnull(patientInfo.iloc[i]['Second_Response_standard'])):
        ResponseType[patientInfo.iloc[i]['Run']]=patientInfo.iloc[i]['Second_Response_standard']
stmp2 = pd.Series(ResponseType)
dict_tmp2 = {'Run':stmp2.index,'Response_type':stmp2.values}
df_tmp2 = pd.DataFrame(dict_tmp2)
patientInfo=pd.merge(patientInfo,df_tmp2,left_on='Run',right_on='Run')
patientInfo=patientInfo.dropna(subset=['Response_type'])
patientInfo.to_excel("./tmp/join_nun_info.xlsx",index=False)
#separate cance type 、anti target and response type
Filter = dict(list(patientInfo.groupby(['Cancer','Anti_target'])))
LungCancer=Filter['lung cancer','anti-PD1'] 
Melanoma_CTLA4=Filter['melanoma','anti-CTLA4']
Melanoma_PD1=Filter['melanoma','anti-PD1']

CR_PR_Lungcancer=LungCancer[['Run']][LungCancer.Response_type=='PR']
PD_Lungcancer=LungCancer[['Run']][LungCancer.Response_type=='PD']
SD_Lungcancer=LungCancer[['Run']][LungCancer.Response_type=='SD']
#Set CR and PR as CR/PR
CR_PR_Lungcancer['Response_type']='CR/PR'   #Response type for box1 is PR
PD_Lungcancer['Response_type']='PD'
SD_Lungcancer['Response_type']='SD'
CR_PR_Lungcancer=CR_PR_Lungcancer.reset_index(drop=True)#重置索引
PD_Lungcancer=PD_Lungcancer.reset_index(drop=True)
SD_Lungcancer=SD_Lungcancer.reset_index(drop=True)

def combine(tmp_dict,tmp_d):#结合同一类sample的信息,返回一个字典
    tmp=tmp_dict[tmp_d.loc[0,"Run"]]
    for s_id in tmp_d["Run"]:
        if (s_id==tmp_d.loc[0,"Run"]):
            continue
        tmp=dict(Counter(tmp)+Counter(tmp_dict[s_id]))
    return (tmp)

    
CR_PR_Hla_dict_lungcancer=combine(hla_dict,CR_PR_Lungcancer)
PD_Hla_dict_dict_lungcancer=combine(hla_dict,PD_Lungcancer)
SD_Hla_dict_dict_lungcancer=combine(hla_dict,SD_Lungcancer)


def get_significant_values(response01_dict,response02_dict): #返回不同类型的有差异的部分
    result={}
    t_list=[]
    for key in response01_dict.keys():
        if (key in response02_dict.keys()):
            if ((response01_dict[key]/response02_dict[key])>1.6 or 
                (response01_dict[key]/response02_dict[key])<0.67):
                t_list.append(response01_dict[key])
                t_list.append(response02_dict[key])
                result[key]=copy.deepcopy(t_list)  #深拷贝，避免重置时带来的错误
        t_list.clear()
    return (result)









