# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:26:29 2017

@author: nmishra function to plot uncertainties and add it to git.
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_smear_signal(smear, int_time_all, quad, color):
    plt.figure()
    for xe, ye in zip(int_time_all, smear):
      plt.scatter([xe] * len(ye), ye, marker='.')  
    
    title = 'Estimated SMEAR Vs. Integration time, '+ quad
    plt.title(title ,fontsize=12,fontweight="bold")
    r" $\mu$" +'secs'
    plt.xlabel('Integration time ('+r"$\mu$" +'secs)')
    plt.ylabel('Smear Signal Estimate (DN)')
    plt.grid(True, linestyle=':')
    plt.ylim(750, 840)
    plt.show()
    #plt.xlim(800, 16000)    
    #plt.savefig(figure_name,dpi=100,bbox_inches="tight")    
    plt.close('all')


def plot_hist_smear_unct(image,quad, COLOR) : 
    
    std_all_cols = 100*np.std(image, axis=0)/np.mean(image, axis=0)
    print(std_all_cols)    
    label = 'Mean = '+ str(round(np.mean(std_all_cols), 2))             
    plt.figure(figsize=(8, 5))
    n, bins, patches = plt.hist(std_all_cols,50, facecolor=COLOR, color='black', label=label, alpha=1)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)    
    plt.xlim(0.1, 0.2)    
    plt.ylim(0, 500)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' SMEAR Signal Uncertainty Estimate (%)' , fontsize=12,
              fontweight="bold")
    plt.title('Histogram of SMEAR Random Uncertainty ('+quad+')')     
    plt.show()
    plt.close('all')
 

def plot_hist_smear(image,quad, COLOR) : 
    
    nx, ny = image.shape
    image = np.reshape(image, (nx*ny,1))
    label = 'Mean = '+ str(round(np.mean(image), 2))             
    plt.figure(figsize=(8, 5))
    n, bins, patches = plt.hist(image,50, facecolor=COLOR, color='black', label=label, alpha=1)
    plt.grid(True, linestyle=':')
    legend = plt.legend(loc='best', ncol=3, shadow=True,
                   prop={'size':10}, numpoints=1)
    legend.get_frame().set_edgecolor('wheat')
    legend.get_frame().set_linewidth(2.0)    
    #plt.xlim(0.1, 0.2)    
    #plt.ylim(0, 700)
    plt.ylabel('Frequency (# of pixels)', fontsize=12,
              fontweight="bold")
    plt.xlabel(' SMEAR Signal Uncertainty Estimate (%)' , fontsize=12,
              fontweight="bold")
    plt.title('Histogram of SMEAR Random Uncertainty ('+quad+')')     
    plt.show()
    plt.close('all')    
    
    
def main():
    
    data_dir = r'C:\Users\nmishra\Workspace\TEMPO\Smear\Using_integ_sweep_data'
    SMEAR_A = np.genfromtxt(data_dir+'/'+'_SMEAR_A.csv', delimiter=',' )
    SMEAR_B = np.genfromtxt(data_dir+'/'+'_SMEAR_B.csv', delimiter=',' )
    SMEAR_C = np.genfromtxt(data_dir+'/'+'_SMEAR_C.csv', delimiter=',' )
    SMEAR_D = np.genfromtxt(data_dir+'/'+'_SMEAR_D.csv', delimiter=',' )
    int_time_all = np.genfromtxt(data_dir+'/'+'all_int_time.csv', delimiter=',' )
    quad = ['Quad A' ,'Quad B', 'Quad C' , 'Quad D']
    colors = ['blue','green','red','orange']
   
#    plot_smear_signal(SMEAR_A, int_time_all, quad[0], colors[0])
#    plot_smear_signal(SMEAR_B, int_time_all, quad[1], colors[1])    
#    plot_smear_signal(SMEAR_C, int_time_all, quad[2], colors[2])    
#    plot_smear_signal(SMEAR_D, int_time_all, quad[3], colors[3])
    
    plot_hist_smear_unct(SMEAR_A, quad[0], colors[0])
    plot_hist_smear_unct(SMEAR_B, quad[1], colors[1])
    plot_hist_smear_unct(SMEAR_C, quad[2], colors[2])
    plot_hist_smear_unct(SMEAR_D, quad[3], colors[3])



if __name__ == "__main__":
    main()
