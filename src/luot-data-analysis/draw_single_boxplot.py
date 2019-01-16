# -*- coding: utf-8 -
'''

'''
import os
import pandas as pd 
import scipy.stats as stats
import xlrd
import matplotlib.pyplot as plt 
import seaborn as sns

pv1=[]
pv2=[]

def get_TILs_Info():#用于获取TILs的信息
    TILs_Info =pd.read_table('../raw_data/TIL_result')
    return TILs_Info

def get_sorted_patients():
    excel=pd.ExcelFile('../raw_data/metadata-verified.xlsx')
    patient_Info=pd.read_excel(excel)
    patient_Info=patient_Info.iloc[:,0:24]#由于最后一行为一些统计的数据所以用切片的方式切掉 
    
    TILs_Info=get_TILs_Info()
    total_sample_name=TILs_Info.iloc[:,0]
    Filter=pd.DataFrame()
    for sample_id in total_sample_name:
        Filter=Filter.append(patient_Info.loc[patient_Info['Run'] == sample_id])#根据Run的列中的ID与TILs中的sample ID匹配出相应的样本信息
    
    Filter.to_excel("../data/First_join_metadata_Tils.xlsx")

    Filter=Filter[Filter['anti-target']=='anti-PD-1']#筛选出anti-target 为anti-PD-1 的个体由于某些缺乏drug信息，所以以anti-target为准。

    Filter = dict(list(Filter.groupby(['Cancer type','response'])))

    get_metastatic_melanoma_R=Filter['metastatic melanoma','R']
    get_metastatic_melanoma_NR=Filter['metastatic melanoma','NR']
    get_metastatic_melanoma_CR=Filter['metastatic melanoma','Complete Response']
    get_metastatic_melanoma_PR=Filter['metastatic melanoma','Partial Response']
    get_metastatic_melanoma_PD=Filter['metastatic melanoma','Progressive Disease']
    get_metastatic_gastric_cancer_CR=Filter['metastatic gastric cancer ','complete response']
    get_metastatic_gastric_cancer_PR=Filter['metastatic gastric cancer ','partial response']
    get_metastatic_gastric_cancer_PD=Filter['metastatic gastric cancer ','progressive disease']
    get_metastatic_gastric_cancer_SD=Filter['metastatic gastric cancer ','stable disease']

    metastatic_melanoma_R=pd.concat([get_metastatic_melanoma_CR,get_metastatic_melanoma_PR,get_metastatic_melanoma_R])
    metastatic_melanoma_NR=pd.concat([get_metastatic_melanoma_PD,get_metastatic_melanoma_NR])
    metastatic_gastric_cancer_R=pd.concat([get_metastatic_gastric_cancer_CR,get_metastatic_gastric_cancer_PR])
    metastatic_gastric_cancer_NR=pd.concat([get_metastatic_gastric_cancer_PD,get_metastatic_gastric_cancer_SD])
    

    return(metastatic_melanoma_R,metastatic_melanoma_NR,metastatic_gastric_cancer_R,metastatic_gastric_cancer_NR)  


def join_TILs_Patients_Info():
    (metastatic_melanoma_R,metastatic_melanoma_NR,metastatic_gastric_cancer_R,metastatic_gastric_cancer_NR)=get_sorted_patients()
    TILs_Info=get_TILs_Info()
    
    Mm_R_Tils=pd.merge(TILs_Info,metastatic_melanoma_R,left_on='sample',right_on='Run')
    Mm_NR_Tils=pd.merge(TILs_Info,metastatic_melanoma_NR,left_on='sample',right_on='Run')
    Mgc_R_Tils=pd.merge(TILs_Info,metastatic_gastric_cancer_R,left_on='sample',right_on='Run')
    Mgc_NR_Tils=pd.merge(TILs_Info,metastatic_gastric_cancer_NR,left_on='sample',right_on='Run')
    return(Mm_R_Tils,Mm_NR_Tils,Mgc_R_Tils,Mgc_NR_Tils)

def write_tilsInfo_to_excel(group1,group2,group3,group4):
    g1=group1.iloc[:,0:25]
    g2=group2.iloc[:,0:25]
    g3=group3.iloc[:,0:25]
    g4=group4.iloc[:,0:25]
    writer=pd.ExcelWriter('../data/all_tils_pInfo.xlsx')
    g1.to_excel(writer,'Sheet1')
    g2.to_excel(writer,'Sheet2')
    g3.to_excel(writer,'Sheet3')
    g4.to_excel(writer,'Sheet4')
    writer.save()

def read_excel(index):
    TILs_R_NR=xlrd.open_workbook('../data/all_tils_pInfo.xlsx')
    cohort_responding=TILs_R_NR.sheet_by_index(0)
    cohort_non_responding=TILs_R_NR.sheet_by_index(1)
    cols_1=cohort_responding.col_values(index)
    cols_2=cohort_non_responding.col_values(index)
    del cols_1[0]
    del cols_2[0]
    return(cols_1,cols_2)

def get_pval(group1,group2):

    s, pVal = stats.ranksums(group1, group2) #  
    return (pVal)

def caculate():
    for i in range(26):
        if i==0 or i==1:
            continue
        (g3,g4)=read_excel(i)
        pv2.append(get_pval(g3,g4))

def DrawSingleBoxplot():
    
    (group1,group2,group3,group4)=join_TILs_Patients_Info()
    g1=group1.iloc[:,1:25]
    g2=group2.iloc[:,1:25]
    g3=group3.iloc[:,1:25]
    g4=group4.iloc[:,1:25]

    c_t=['CD4_naive', 'CD8_naive', 'Cytotoxic', 'Exhausted', 'Tr1',
       'nTreg', 'iTreg', 'Th1', 'Th2', 'Th17', 'Tfh', 'Central_memory',
       'Effector_memory', 'NKT', 'MAIT', 'DC', 'Bcell', 'Monocyte',
       'Macrophage', 'NK', 'Neutrophil', 'Gamma_delta', 'CD4_T', 'CD8_T']

    for i in range(24):

        gr_1=pd.DataFrame()
        gr_2=pd.DataFrame()
        gr_3=pd.DataFrame()
        gr_4=pd.DataFrame()

        tmp1=g1.iloc[:,i:i+1]
        tmp1.columns=['amount']
        tmp1['response_type']='R'
        tmp1['cell_type']=c_t[i]
        gr_1=gr_1.append(tmp1)
        

        tmp2=g2.iloc[:,i:i+1]
        tmp2.columns=['amount']
        tmp2['response_type']='NR'
        tmp2['cell_type']=c_t[i]
        gr_2=gr_2.append(tmp2)
        
        tmp3=g3.iloc[:,i:i+1]
        tmp3.columns=['amount']
        tmp3['response_type']='R'
        tmp3['cell_type']=c_t[i]
        gr_3=gr_3.append(tmp3)

        tmp4=g4.iloc[:,i:i+1]
        tmp4.columns=['amount']
        tmp4['response_type']='NR'
        tmp4['cell_type']=c_t[i]
        gr_4=gr_4.append(tmp4)

        
    
        mm_join=pd.concat([gr_1,gr_2])
        mgc_join=pd.concat([gr_3,gr_4])
       
        mm_join.to_excel('../data/mm_boxplot.xlsx')
        mgc_join.to_excel('../data/mgc_boxplot.xlsx')

        plt.figure(figsize=(7, 10))
        plt.title("p="+str('%.3f' % pv2[i]))  #pv2[i]
        
        ax = sns.boxplot(x="cell_type", y="amount", hue="response_type",data=mm_join, palette="Set3",fliersize=0) #data=mm_join
        plt.savefig("../boxplot/mm/"+str(i+1)+"_"+c_t[i]+".pdf") #/mm mgc

def main():
    caculate()
    DrawSingleBoxplot()
    
    

if __name__ == '__main__':
    main()
