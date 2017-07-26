# microbiome-analysis
#how can we understand mircobiome
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
import skbio
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes
from numpy import std, mean, sqrt
%matplotlib inline


#Alpha diversity analysis
def Alpha_analysis(alpha_df): 
    #input data frame contains columns of shannon index, following by feature of each sample. Row names are sample ID.
    #create a feature list
    del alpha_df['Observed OTU']
    features_df=alpha_df
    features_df=features_df.drop(features_df.columns[0], axis=1)
    n = features_df.shape[1]
    features_list = [[] for _ in range(n)]
    for i in range(n):
        features_list[i]=(list(set(features_df.ix[:,i].values)))
        
    #split data into different features
    lists = [[] for _ in range(n*2)]
    p_list = [[] for _ in range(n)]
    for i in range (len(features_list)):
        for j in range (alpha_df.shape[0]):
            if alpha_df.ix[j][i+1] == features_list[i][0]:
                lists[i*2].append(alpha_df.ix[j][0])
            else:
                lists[i*2+1].append(alpha_df.ix[j][0])
    #find p value
    for i in range (n):
        z,p=scipy.stats.ranksums(lists[i*2],lists[i*2+1])
        p_list[i].append(p)
    p_df = pd.DataFrame(data=p_list, index=features_df.columns.values, columns=['p'])
    p_df.to_csv('/home/gene/sunzheng/pm_duplicate/p_alpha.csv', sep='\t')
    #plot alpha diversity
    for i in range (n):
        feature_1=pd.DataFrame(data=lists[i*2], index=None, columns=[features_list[i][0]])
        feature_2=pd.DataFrame(data=lists[i*2+1], index=None, columns=[features_list[i][1]])
        data_df=pd.concat([feature_1, feature_2], axis=1)
        figure(i+1)
        sns.set(font_scale=1.5)
        ax1=sns.boxplot(data_df)
        plt.xlabel(alpha_df.columns.values[i+1], fontsize=20)
        plt.ylabel('Shannon', fontsize=20)
        plt.savefig('/home/gene/sunzheng/pm_duplicate/'+'alpha' + alpha_df.columns.values[i+1] +'.pdf')

#Beta diversity analysis
def beta_analysis(distance_df,feature_df):
    #match features with each distance data
    n=feature_df.shape[1]
    data=[]
    lists_featureid = [[] for _ in range(n*2)]
    for i in range (distance_df.shape[0]-1):
        for j in range (i+1,distance_df.shape[1]):
            data.append(distance_df.ix[i,j])
            for z in range(n):
                lists_featureid[z*2].append(feature_df.ix[i][z])
                lists_featureid[z*2+1].append(feature_df.ix[j][z])
    
    #create a feature list 
    feature_list = [[] for _ in range(n)]
    for i in range(n):
        feature_list[i]=(list(set(feature_df.ix[:,i].values)))
    
    #split data for each feature in within or between
    data_featuregroup = [[] for _ in range(n*2)]
    for i in range(n):
        for j in range (len(data)):
            if lists_featureid[i*2][j]==lists_featureid[i*2+1][j]:
                data_featuregroup[i*2].append(data[j])
            else:
                data_featuregroup[i*2+1].append(data[j])
    
    #effect size
    def cohen_d(x,y):
        nx = len(x)
        ny = len(y)
        dof = nx + ny - 2
        return (mean(x) - mean(y)) / sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)
    
    #plot data
    # function for setting the colors of the box plots pairs
    def setBoxColors(bp):
        setp(bp['boxes'][0], color='blue')
        setp(bp['caps'][0], color='blue')
        setp(bp['caps'][1], color='blue')
        setp(bp['whiskers'][0], color='blue')
        setp(bp['whiskers'][1], color='blue')
        setp(bp['fliers'][0], markeredgecolor='blue')
        setp(bp['medians'][0], color='blue')

        setp(bp['boxes'][1], color='red')
        setp(bp['caps'][2], color='red')
        setp(bp['caps'][3], color='red')
        setp(bp['whiskers'][2], color='red')
        setp(bp['whiskers'][3], color='red')
        setp(bp['fliers'][1], markeredgecolor='red')
        setp(bp['medians'][1], color='red')
    
    fig = figure()
    ax = axes()
    hold(True)
    sns.set(font_scale=1.5)
    d=[]
    for i in range(n):
        data_plot=[data_featuregroup[i*2],data_featuregroup[i*2+1]]
        bp = boxplot(data_plot, positions = [1+3*i, 2+3*i], widths = 0.6)
        setBoxColors(bp)
        effect_size=cohen_d(data_featuregroup[i*2],data_featuregroup[i*2+1])
        d.append(effect_size)

    def floatrange(start,stop,steps):
        return [start+float(i)*(stop-start)/(float(steps)-1) for i in range(steps)]
    # set axes limits and labels
    xlim(0,n*3)
    ylim(0,max(data)+0.05)
    ax.set_xticklabels(feature_df.columns.values) 
    ax.set_xticks(floatrange(1.5,n*3-1.5,3))
    ax.set_ylabel('Distance') 


    # draw temporary red and blue lines and use them to create a legend
    hB, = plot([1,1],'b-')
    hR, = plot([1,1],'r-')
    legend((hB, hR),('Within', 'Between'))
    hB.set_visible(False)
    hR.set_visible(False)
    
    #beta diversity p value from permanova
    #two inputs, distance data and grouping 
    distance_data = distance_df.values
    distance_ids = distance_df.index.values
    distance_matrix = skbio.stats.distance.DistanceMatrix(distance_data, distance_ids)
    p_list=[]
    for i in range(n):
        grouping = feature_df.ix[:,i].values
        p_beta=skbio.stats.distance.permanova(distance_matrix, grouping, permutations=999)['p-value']
        p_list.append(p_beta)
    
    #save file 
    #plt.savefig('/home/gene/sunzheng/pm_duplicate/beta_pt.pdf')
    p_df = pd.DataFrame(data=p_list, index=feature_df.columns.values, columns=['p'])
    e_df = pd.DataFrame(data=d, index=feature_df.columns.values, columns=['Cohens_d'])
    result_df=pd.concat([p_df, e_df], axis=1)
    #result_df.to_csv('/home/gene/sunzheng/pm_duplicate/beta_data.csv', sep='\t')


#Bio-marker finding
def biomarker(abundance_df,feature_df):
    #input is the abundance table. The last column contain the feature infor
    feature_list = list(set(feature_df.ix[:,-1].values))
    count_matrix=np.zeros((abundance_df.shape[1], 1))
    for i in range (abundance_df.shape[1]):
        feature = [[] for _ in range(2)]
        for j in range (abundance_df.shape[0]):
            if feature_df.ix[j][-1] == feature_list[0]:
                feature[0].append(abundance_df.ix[j][i])
            else:
                feature[1].append(abundance_df.ix[j][i])
            z_stat, p_val=scipy.stats.ranksums(feature[0],feature[1])
            count_matrix[i]=p_val
    index=abundance_df.columns.values
    sort_df=pd.DataFrame(data=count_matrix, index=index, columns=['P_val'])
    sort_df=sort_df[sort_df['P_val']<=0.01].sort_values('P_val')
    sort_df.to_csv('/home/gene/sunzheng/pm_duplicate/biomarker.csv', sep='\t')
    
    #find genus with higher relative abundance 
    threshold_id=abundance_df[abundance_df>0.02]
    threshold_id=threshold_id.dropna(axis=1, thresh=abundance_df.shape[0]/2).columns.values
    
    #plot marker with high relative abundance
    for i in range(len(threshold_id)):
        f1=[]
        f2=[]
        if threshold_id[i] in sort_df.index.values:
            ab_data=abundance_df[threshold_id[i]]
            for j in range(len(ab_data)):
                if feature_df.ix[:,-1][j]==feature_list[0]:
                    f1.append(ab_data[j])
                else:
                    f2.append(ab_data[j])
            f1=pd.DataFrame(data=f1, index=None, columns=[feature_list[0]])
            f2=pd.DataFrame(data=f2, index=None, columns=[feature_list[1]])
            data_df=pd.concat([f1, f2], axis=1)
            figure(i+1)
            sns.set(font_scale=1.5)
            ax1=sns.boxplot(data_df,vert=False)
            plt.xlabel('Relative Abundance', fontsize=20)
            plt.ylabel(feature_df.columns.values[0], fontsize=20)
            ax1.set_title(threshold_id[i],fontsize=20)
            plt.savefig('/home/gene/sunzheng/pm_duplicate/'+'biomarker' + threshold_id[i] +'.pdf')    
