# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sst
import math
import time

time_start = time.time()
c5_list = [330]
#color_list = ['yellow', 'red', 'blue']
color_list = ['red', 'red', 'blue']
label_list =['Chimera', 'Coherent', 'Metastable']

chimera_id = np.load("data_results_sigma_001_1/WCOs_index_chimera_c5_" + str(c5_list[0])+".npy")

n_trial = 10

x0 = 0
y0 = 0
rho_sst = []
C_sst = []

plt.figure(figsize=[6.4*1.4, 4.8*1.4])
for i in range(len(c5_list)):
    c5 = c5_list[i]
    C_arr = np.load("data_results_sigma_001_1/WCOs_chimeta_id_"+str(c5)+"_aal2_N94_0.npy")
    C_arr = C_arr[chimera_id[i]]
    x0 += len(C_arr)
    C_sst += C_arr.tolist()
    results_rho = np.load("data_results_sigma_001_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_0.npy")
    results_rho = results_rho[chimera_id[i]]
    y0 += len(results_rho)
    rho_sst += results_rho.tolist()
    plt.plot(results_rho, C_arr, '.', color=color_list[i], markersize=10,label=label_list[i] ,alpha=0.8)
    for j in np.arange(1, n_trial):
        C_arr = np.load("data_results_sigma_001_1/WCOs_chimeta_id_" + str(c5) + "_aal2_N94_"+str(j)+".npy")
        C_arr = C_arr[chimera_id[j]]
        x0 += len(C_arr)
        C_sst += C_arr.tolist()
        results_rho = np.load("data_results_sigma_001_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_"+str(j)+".npy")
        results_rho = results_rho[chimera_id[j]]
        y0 += len(results_rho)
        rho_sst += results_rho.tolist()
        plt.plot(results_rho, C_arr, '.', color=color_list[i], markersize=10, alpha=0.8)
#plt.title(r"$c_5=$"+str(c5_list), fontsize=15)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r"$r_N$", fontsize=20)
plt.ylabel(r"$\Gamma$", fontsize=20)
#plt.legend(fontsize=20)
plt.savefig("Figures_sigma_001_1/chimera_id_vs_global_synch_1.pdf")

print(sst.pearsonr(rho_sst, C_sst))

plt.show()


plt.figure(figsize=[6.4*1.1, 4.8*1.1])
i=2
c5 = 200
lambda_arr = np.load("data_results_sigma_001_1/WCOs_metastability_id_"+str(c5)+"_aal2_N94_0.npy")
results_rho = np.load("data_results_sigma_001_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_0.npy")
plt.plot(results_rho, lambda_arr, '.', color=color_list[i], markersize=10,label=label_list[i] ,alpha=0.8)
for j in np.arange(1, n_trial):
    lambda_arr = np.load("data_results_sigma_001_1/WCOs_metastability_id_" + str(c5) + "_aal2_N94_"+str(j)+".npy")
    results_rho = np.load("data_results_sigma_001_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_"+str(j)+".npy")
    plt.plot(results_rho, lambda_arr, '.', color=color_list[i], markersize=10, alpha=0.8)
#plt.title(r"$c_5=$"+str(c5), fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r"Global synchrony", fontsize=15)
plt.ylabel(r"Metastabiliy idex", fontsize=15)
#plt.legend(fontsize=15)
plt.savefig("Figures_sigma_001_1/meta_id_vs_global_synch_1.pdf")
plt.show()

time_terminal = time.time()
print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")