# -*- coding: utf-8 -*-
from __future__ import division
import os
import re
import glob
import copy
import xlwt
import numpy as np
import pandas as pd 
from collections import Counter
from ShowProcess import ShowProcess




patientInfo=pd.read_csv('clinical_info',sep='\t',skipinitialspace=True )
patientInfo=patientInfo.loc[patientInfo['Biopsy_Time']=='pre-treatment']
sample_list=pd.read_csv('./neoantigen_filter_SB_result/sample_list')   ###################记得修改
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




h_n_p_list=[]
h_n_p_dict={}

filelist=(glob.glob(r'./neoantigen_filter_SB_result/*.filter'))      ###################记得修改

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
            s=''
            line = re.sub(' +', '#', line)  #解决分隔符不统一的问题，将出现空格的地方替换为 “#”
            line=line.split("#")
            s+=line[19]         #添加基因信息
            s+=","              
            s+=line[20]         #添加新抗原信息
            s+="," 
            s+=line[28]         #添加新的蛋白质的信息
            h_n_p_list.append(s)

        if  sample_id in h_n_p_dict.keys():
            temp1 = all_np(h_n_p_list)
            temp2 = h_n_p_dict[sample_id]    #这是一个字典  充分利用python的浅拷贝原理改变字典的值
            for k in np.unique(h_n_p_list):
                if (k not in temp2.keys()):
                    temp2[k] = temp1[k]
                temp2[k] += temp1[k]
            h_n_p_list.clear()
        else:
            h_n_p_dict[sample_id]=all_np(h_n_p_list)
            h_n_p_list.clear()

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

count_LungCancer=len(LungCancer)   
count_Melanoma_CTLA4=len(Melanoma_CTLA4)
count_Melanoma_PD1=len(Melanoma_PD1)


CR_PR_Lungcancer=LungCancer[['Run']][LungCancer.Response_type=='PR']
PD_Lungcancer=LungCancer[['Run']][LungCancer.Response_type=='PD']
SD_Lungcancer=LungCancer[['Run']][LungCancer.Response_type=='SD']
LTB_Melanoma_CTLA4=Melanoma_CTLA4[['Run']][Melanoma_CTLA4.Response_type=='long-term-benefit']
MONB_Melanoma_CTLA4=Melanoma_CTLA4[['Run']][Melanoma_CTLA4.Response_type=='minimal-or-no-benefit']
CR_PR_Melanoma_PD1=Melanoma_PD1[['Run']][(Melanoma_PD1.Response_type=='CR') |(Melanoma_PD1.Response_type=='PR') ]
PD_Melanoma_PD1=Melanoma_PD1[['Run']][Melanoma_PD1.Response_type=='PD' ]


#Set CR and PR as CR/PR
CR_PR_Lungcancer['Response_type']='CR/PR'   #Response type for box1 is PR
PD_Lungcancer['Response_type']='PD'
SD_Lungcancer['Response_type']='SD'
LTB_Melanoma_CTLA4['Response_type']='long-term-benefit'
MONB_Melanoma_CTLA4['Response_type']='minimal-or-no-benefit'
CR_PR_Melanoma_PD1['Response_type']='CR/PR'
PD_Melanoma_PD1['Response_type']='PD'


def combine(tmp_dict,tmp_d):#结合同一类sample的信息,返回一个字典。
    tmp=tmp_dict[tmp_d.iloc[0,0]]
    for tmpKey in tmp.keys():
        tmp[tmpKey] = [tmp[tmpKey] ,[tmp_d.iloc[0,0]]]
    for s_id in tmp_d["Run"]:
        if (s_id==tmp_d.iloc[0,0]):
            continue
        if (s_id in tmp_dict.keys()): #某些文件为空，所以部分Sample id 会有可能找不到
            a = Counter(tmp_dict[s_id])    #将重复的累计起来
            flag = 0
            for  aKey in a.keys():
                if (aKey in tmp.keys()):
                    flag += 1
                    tmp[aKey][0] += a[aKey]
                    if (flag == 1):
                        tmp[aKey][1].append(s_id)
                else:
                    tmp[aKey] = [0,[]]
                    tmp[aKey][0] = a[aKey]
                    tmp[aKey][1].append(s_id)
    return (tmp)

CR_PR_dict_lungcancer=combine(h_n_p_dict,CR_PR_Lungcancer)
PD_dict_lungcancer=combine(h_n_p_dict,PD_Lungcancer)
SD_dict_lungcancer=combine(h_n_p_dict,SD_Lungcancer)

LTB_dict_Melanoma_CTLA4=combine(h_n_p_dict,LTB_Melanoma_CTLA4)
MONB_dict_Melanoma_CTLA4=combine(h_n_p_dict,MONB_Melanoma_CTLA4)

CR_PR_dict_Melanoma_PD1=combine(h_n_p_dict,CR_PR_Melanoma_PD1)
PD_dict_Melanoma_PD1=combine(h_n_p_dict,PD_Melanoma_PD1)




def get_significant_values(response01_dict,response02_dict): #返回不同类型的有差异的部分
    result={}
    t_list=[]
    ratio_dict={}
    for key in response01_dict.keys():
        if (key in response02_dict.keys()):


            if ((response01_dict[key][0]/response02_dict[key][0])>1.5 or (response01_dict[key][0]/response02_dict[key][0])<0.67):
                t_list.append(response01_dict[key][0])
                t_list.append(response02_dict[key][0])
                result[key]=copy.deepcopy(t_list)  #深拷贝，避免重置时带来的错误
                ratio_dict[key]=[]
                ratio_dict[key].append(len(response01_dict[key][1])+len(response02_dict[key][1]))
                ratio_dict[key].append(len(response01_dict[key][1]))
                ratio_dict[key].append(len(response02_dict[key][1]))
                t_list.clear()
    

    return (result,ratio_dict)


lungcancer_cr_pr_pd,count_lungcancer_cr_pr_pd=get_significant_values(CR_PR_dict_lungcancer,PD_dict_lungcancer)


lungcancer_cr_pr_sd,count_lungcancer_cr_pr_sd=get_significant_values(CR_PR_dict_lungcancer,SD_dict_lungcancer)

Melanoma_CTLA4_ltb_monb,count_Melanoma_CTLA4_ltb_monb=get_significant_values(LTB_dict_Melanoma_CTLA4,MONB_dict_Melanoma_CTLA4)

Melanoma_PD1_cr_pr_pd,count_Melanoma_PD1_cr_pr_pd=get_significant_values(CR_PR_dict_Melanoma_PD1,PD_dict_Melanoma_PD1)


def write_to_excel(t_dict,name,r1,r2,total,count):  #将结果写入Excel文档
    wbk = xlwt.Workbook()
    sheet = wbk.add_sheet('sheet 1')
    sheet.write(0,1,r1)
    sheet.write(0,2,r2)
    sheet.write(0,3,"Total ratio")
    sheet.write(0,4,r1+" ratio")
    sheet.write(0,5,r2+" ratio")
    for ca,i in zip(t_dict.keys(),range(100)):
        sheet.write(i+1,0,ca)
        sheet.write(i+1,1,t_dict[ca][0])
        sheet.write(i+1,2,t_dict[ca][1])
        print(count[ca])
        sheet.write(i+1,3,count[ca][0]/total)
        sheet.write(i+1,4,count[ca][1]/total)
        sheet.write(i+1,5,count[ca][2]/total)
    wbk.save("./result/"+name+'.xls')
    print(name+" write to excel!")


write_to_excel(lungcancer_cr_pr_pd,"lungcancer_cr_pr_pd","CR/PR","PD",count_LungCancer,count_lungcancer_cr_pr_pd)
write_to_excel(lungcancer_cr_pr_sd,"lungcancer_cr_pr_sd","CR/PR","SD",count_LungCancer,count_lungcancer_cr_pr_sd)

write_to_excel(Melanoma_CTLA4_ltb_monb,"Melanoma_CTLA4_ltb_monb","long-term-benefit","minimal-or-no-benefit",count_Melanoma_CTLA4,count_Melanoma_CTLA4_ltb_monb)

write_to_excel(Melanoma_PD1_cr_pr_pd,"Melanoma_PD1_cr_pr_pd","CR/PR","PD",count_Melanoma_PD1,count_Melanoma_PD1_cr_pr_pd)
