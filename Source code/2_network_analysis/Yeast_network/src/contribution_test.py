#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 15:02:19 2025

@author: caio
"""

from metrics_lib import Get_metrics_data, adapt_yeast
import numpy as np
from sklearn.datasets import load_breast_cancer, load_iris
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import os

def Feature_Contribution_Matrix(X):
  
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    
    pca = PCA().fit(X_scaled) 
    X_new = PCA().fit_transform(X) # project the original data into the PCA space
    

    dot_matrix = np.dot(X.T, X_new)
   
    dot_matrix_norm = dot_matrix.copy()
    
    for j in range(pca.n_components_):
        col = dot_matrix_norm[:,j:j+1]   
        dot_matrix_norm[:,j:j+1] = (col - col.mean()) / col.std(ddof = 1)
    
    exp_var = pca.explained_variance_ratio_
    norm_corr = np.abs(dot_matrix_norm)*exp_var
    
    return dot_matrix_norm, norm_corr, pca


def Plot_Heatmap(heatmap, ax,
                 xlabels, ylabels,
                 title = '', cmap = 'Purples', fs = 25):


    imshow = ax.imshow(heatmap,
               cmap = cmap)
    ax.set_title(title, size = fs+10)
    
    ax.set_xticks(np.arange(len(xlabels)),
                  xlabels, size = fs-5,
                  rotation = 45, rotation_mode = 'anchor',
                  horizontalalignment='right', verticalalignment='top')
    
    ax.set_yticks(np.arange(len(ylabels)),
                  ylabels, size = fs-5)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mappable = imshow,
                 cax = cax).ax.tick_params(labelsize=fs-5)



def prep_plot(num_rows, num_cols, scale_fig=8, ratio = [1,1]):
    
    """
    ratio = [rx, ry]
        width  ->  rx*scale_fig*num_cols
        height ->  ry*scale_fig*num_rows
    """
    
    

    fig = plt.figure(figsize = (ratio[0]*scale_fig*num_cols,  #width
                                ratio[1]*scale_fig*num_rows)) #height
           

    ax_list = []
    for i in range(num_cols*num_rows):
        ax_list.append(fig.add_subplot(num_rows, num_cols, i+1))

    return fig


def Print_Variance_PCA(pca):
    # Get cumulative explained variance for each dimension
    pca_evr = pca.explained_variance_ratio_
    cumsum_ = np.cumsum(pca_evr)
    
    # Get dimensions where var >= 95% and values for variance at 2D, 3D
    dim_95 = np.argmax(cumsum_ >= 0.95) + 1
        
    dims_ = pca.n_components_

    
    # Print report
    text = ''
    text+= f"You can reduce from {dims_} to {dim_95} dimensions while retaining 95% of variance."
    text+='\n'
    for i in range(dims_):
        i_dim = np.round(cumsum_[i], decimals=3)*100 
        text+=f"  {i+1} principal components explain {i_dim:.2f}% of variance."
        text+='\n'    
    print(text)
    print()
    return text



def Four_Plots(df_data, fs = 15):
    
    fig = prep_plot(2,2,ratio = [1,1], scale_fig = 20 )
    axes = fig.get_axes()
    
    features = list(df_data.columns)
    data = np.array(df_data)
    #Plotting the correlation matrix of the features
    correlation_matrix = df_data.corr()
    ax0 = axes[0]
    Plot_Heatmap(correlation_matrix, ax0,
                 features, features, 
                 title = 'Correlation Matrix',
                 cmap='Reds', fs = fs)
    
    ax0.set_xticks(np.arange(len(features)), features,
                   rotation = 'vertical', size = fs)
    
    
    dot_matrix_norm, norm_corr, pca = Feature_Contribution_Matrix(data)
    text = Print_Variance_PCA(pca)
    
    X_pca = pca.fit_transform(data)
    ax1 = axes[1]
    ax1.scatter(X_pca[:,0],
                X_pca[:,1],
                s = 70,
                c = 'darkgreen')
    ax1.set_xlabel('PC1', fontsize = fs)
    ax1.set_ylabel('PC2', fontsize = fs)
    ax1.set_title('Principal Components', fontsize = fs+10)
    ax1.tick_params(labelsize = fs)
    
        
    
    PCA_labels = [f'PC{i+1}' for i in range(pca.n_components_)]
    ax2 = axes[2]
    Plot_Heatmap(dot_matrix_norm, ax2,
                 PCA_labels, features, 
                 title = 'Feature Contribution',
                 cmap = 'PiYG', fs = fs)
    
    
    
    ax3 = axes[3]
    Plot_Heatmap(norm_corr, ax3,
                 PCA_labels, features, 
                 title = '|Feature Contribution|*Explained Variance ',
                 cmap='PuBuGn', fs = fs)
    
    
    plt.tight_layout()
    
    return fig, text
    
    
    
def PCA_plot(df_data, ax = plt.gca(),
             title = '', cmap = 'PuBuGn',
             do_print = False, fs = 15):
    
    
    data = np.array(df_data)
    features = list(df_data.columns)

    dot_matrix_norm, norm_corr, pca = Feature_Contribution_Matrix(data)

    PCA_labels = [f'PC{i+1}' for i in range(pca.n_components_)]

    Plot_Heatmap(norm_corr, ax,
                 PCA_labels, features, 
                 title = title,
                 cmap = cmap, fs = fs)
    
    results = []
    for i in range(len(features)):
        if do_print:
            print(features[i], norm_corr[i].sum())
        results.append([features[i], norm_corr[i].sum()])
        
    results = pd.DataFrame(results,
                           columns = ['Metric', 'Contribution(sum of PCA component)'])
    
    results.set_index('Metric', inplace = True)

                           
    return results
    
    
    
    
if __name__ == "__main__":

    # graph = 'trpA'
    # metrics = Get_metrics_data(graph)


    # PCA_plot(metrics)


    # fig4P, text = Four_Plots(metrics, fs = 25)
    
    

    # data_folder = '../data/subnets_table'
    # subnets = os.listdir(data_folder)
    # graph = subnets[0]


    # metrics, gabarito = adapt_yeast(graph, data_folder)
    # metrics.pop('gabarito')
    # results = PCA_plot(metrics)


    # # fig4P, text = Four_Plots(metrics, fs = 25)
    
    
    
    
    
    data_folder = '../metanalysis/chimeras'
    subnets = os.listdir(data_folder)
    graph = 'CL_137_18.12_combined_score.csv'

    metrics = Get_metrics_data(graph, data_folder)

    results = PCA_plot(metrics)

    