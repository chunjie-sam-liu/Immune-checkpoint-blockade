import os
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

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
Filter = dict(list(patientInfo.groupby(['Cancer','Anti_target'])))
LungCancer=Filter['lung cancer','anti-PD1'] 
Melanoma_CTLA4=Filter['melanoma','anti-CTLA4']
Melanoma_PD1=Filter['melanoma','anti-PD1']

def draw_boxplot_lungcancer():
    box1=(LungCancer[['MutNumber']][LungCancer.Response_type=='PR']).iloc[:,0].values
    box2=(LungCancer[['MutNumber']][LungCancer.Response_type=='PD']).iloc[:,0].values
    box3=(LungCancer[['MutNumber']][LungCancer.Response_type=='SD']).iloc[:,0].values

    s, pVal = stats.mannwhitneyu(box3, box1) #  wilcox秩序和检验ranksums Mann-Whitney U检验 mannwhitneyu                             
    d=[box1,box2,box3]
    fig, axes = plt.subplots(figsize=(4, 7))
    fig.subplots_adjust(left=0.25,right=0.90,bottom=0.11
            ,top=0.88,wspace=0.2,hspace=0.2 )
    bplot=axes.boxplot(d,labels=['PR','PD','SD'],patch_artist=True)
    axes.set_title('lung cancer')
    axes.set_xlabel('Response type')
    axes.set_ylabel('Mutation number')

    colors = ['lightgreen', 'yellow', 'pink']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    plt.savefig("lung_cancer.pdf")


def draw_boxplot_Melanoma_CTLA4():
    box1=Melanoma_CTLA4[['MutNumber']][Melanoma_CTLA4.Response_type=='long-term-benefit'].iloc[:,0].values
    box2=Melanoma_CTLA4[['MutNumber']][Melanoma_CTLA4.Response_type=='minimal-or-no-benefit'].iloc[:,0].values
    d=[box1,box2]

    s, pVal1 = stats.mannwhitneyu(box2, box1)
    print(pVal1)
    s, pVal2 = stats.ranksums(box2, box1)
    print(pVal2)

    fig, axes = plt.subplots(figsize=(4, 7))
    fig.subplots_adjust(left=0.22,right=0.90,bottom=0.25
                ,top=0.95,wspace=0.2,hspace=0.2 )
    bplot=axes.boxplot(d,labels=['long-term-benefit','minimal-or-no-benefit'],patch_artist=True)
    axes.set_title('Melanoma_CTLA4')
    axes.set_xlabel('Response type')
    axes.set_ylabel('Mutation number')
    plt.xticks(rotation=60)
    colors = ['lightgreen', 'pink']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    plt.savefig("Melanoma_CTLA4.pdf")

def draw_boxplot_Melanoma_PD1():
    box1=Melanoma_PD1[['MutNumber']][(Melanoma_PD1.Response_type=='CR') |(Melanoma_PD1.Response_type=='PR') ].iloc[:,0].values
    box2=Melanoma_PD1[['MutNumber']][Melanoma_PD1.Response_type=='PD' ].iloc[:,0].values

    s, pVal1 = stats.mannwhitneyu(box2, box1)
    print(pVal1)
    s, pVal2 = stats.ranksums(box2, box1)
    print(pVal2)
    
    d=[box1,box2]
    fig, axes = plt.subplots(figsize=(4, 7))
    fig.subplots_adjust(left=0.21,right=0.90,bottom=0.08
                ,top=0.88,wspace=0.2,hspace=0.2 )
    bplot=axes.boxplot(d,labels=['CR/PR','PD'],patch_artist=True)
    axes.set_title('Melanoma_PD1')
    axes.set_xlabel('Response type')
    axes.set_ylabel('Mutation number')
    colors = ['lightgreen', 'pink']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
    plt.savefig("Melanoma_PD1.pdf")

def main():
    pass

if __name__ == '__main__':
    main()




