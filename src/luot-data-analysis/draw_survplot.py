# -*- coding: utf-8 -
'''
    为画生存图做准备
'''
import os
import pandas as pd 
import scipy.stats as stats
import xlrd
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np

c_t=['CD4_naive', 'CD8_naive', 'Cytotoxic', 'Exhausted', 'Tr1',
       'nTreg', 'iTreg', 'Th1', 'Th2', 'Th17', 'Tfh', 'Central_memory',
       'Effector_memory', 'NKT', 'MAIT', 'DC', 'Bcell', 'Monocyte',
       'Macrophage', 'NK', 'Neutrophil', 'Gamma_delta', 'CD4_T', 'CD8_T']

def GetTilsInfo():#用于获取TILs的信息
    TILs_Info =pd.read_table('../raw_data/TIL_result')
    return TILs_Info

def GetSortedPatients():
    excel=pd.ExcelFile('../raw_data/metadata-verified.xlsx')
    patient_Info=pd.read_excel(excel)
    patient_Info=patient_Info.iloc[:,0:24]#由于最后一行为一些统计的数据所以用切片的方式切掉 
    
    TILs_Info=GetTilsInfo()
    total_sample_name=TILs_Info.iloc[:,0]
    Filter=pd.DataFrame()
    for sample_id in total_sample_name:
        Filter=Filter.append(patient_Info.loc[patient_Info['Run'] == sample_id])#根据Run的列中的ID与TILs中的sample ID匹配出相应的样本信息
    Filter = dict(list(Filter.groupby(['Cancer type'])))
    Filter=Filter['metastatic melanoma']
    Filter=Filter[['Cancer type','Run','Survival time(day)','Survival status']] #选取相应的有用的数据
    Filter=Filter.dropna(axis=0,subset=['Survival status'])#过滤掉没有Survival status信息的样本
    Filter=pd.merge(TILs_Info,Filter,left_on='sample',right_on='Run')
    return Filter


def SortBottomTopTils():
    tils_surv=GetSortedPatients()
    j=1
    for i in c_t:
        tmp=tils_surv[[i,'Cancer type','Run','Survival time(day)','Survival status']]
        tmp=tmp.sort_index(axis = 0,ascending = True,by =i)
        bottom=tmp[0:9]
        bottom['type']='bottom third'
        top=tmp[19:]
        top['type']='top third'
        bottom_top_sorted=pd.concat([bottom,top])
        bottom_top_sorted = bottom_top_sorted.rename(columns={'Survival time(day)':'Survival_time','Survival status':'Survival_status'})
        bottom_top_sorted.to_csv('../../survplot/third/'+str(j)+"_"+i+".csv",index=False)
        print("Done!")
        j=j+1


def main():
    SortBottomTopTils()

if __name__ == '__main__':
    main()