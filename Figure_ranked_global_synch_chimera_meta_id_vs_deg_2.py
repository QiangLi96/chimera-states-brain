# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import optimize
import scipy.stats as sst
import time

time_start = time.time()
A_mat = np.load("data3T_aal2_A_count_N94_1.npy")
deg_mat = np.sum(A_mat, axis=1)
N = 94
c5 = 330
n_trial = 10

"""
plt.figure(figsize=[6.4 * 1.6, 4.8 * 1.4])
rho_list = []
deg_list = []
for i in range(n_trial):
    #C_arr = np.load("WCOs_chimeta_id_"+str(c5)+"_aal2_N94_1.npy")
    results_rho = np.load("data_results_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_"+str(i)+".npy")
    rho_list.append(results_rho)
    deg_list.append(deg_mat)
    plt.plot(deg_mat, results_rho, '.', color='gray', markersize=15, alpha=0.8)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r"Ranked weighted degree", fontsize=15)
    plt.ylabel(r"Ranked global synchrony", fontsize=15)

rho_list = np.array(rho_list)
rho_list = rho_list.flatten('C')
deg_list = np.array(deg_list)
deg_list = deg_list.flatten('C')
print(sst.pearsonr(deg_list,rho_list))

plt.figure(figsize=[6.4 * 1.6, 4.8 * 1.4])
C_list = []
for i in range(n_trial):
    C_arr = np.load("data_results_1/WCOs_chimeta_id_" + str(c5) + "_aal2_N94_"+str(i)+".npy")
    C_list.append(C_arr)
    plt.plot(deg_mat, C_arr, '.', color='gray', markersize=15, alpha=0.8)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r"Ranked weighted degree", fontsize=15)
    plt.ylabel(r"Ranked chimera index", fontsize=15)
C_list = np.array(C_list)
C_list = C_list.flatten('C')
print(sst.pearsonr(deg_list,C_list))

"""

chimera_id = np.load("data_results_sigma_001_1/WCOs_index_chimera_c5_" + str(c5)+".npy")




#直线方程函数
def f_1(x, a, b):
    return a*x + b


deg_rank = []
rho_rank = []
x0 = 0
y0 = 0

plt.figure(figsize=[6.4 * 1.2, 4.8 * 1.4])
for i in range(n_trial):
    # C_arr = np.load("WCOs_chimeta_id_"+str(c5)+"_aal2_N94_1.npy")
    results_rho = np.load("data_results_sigma_001_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_"+str(i)+".npy")
    rho_tmp = results_rho[chimera_id[i]]
    deg_tmp = deg_mat[chimera_id[i]]
    rho_tmp = list(rho_tmp)
    deg_tmp = list(deg_tmp)
    x0 += len(rho_tmp)
    y0 += len(deg_tmp)

    id_deg = np.argsort(np.argsort(deg_tmp))
    deg_rank.append(id_deg)
    id_rho = np.argsort(np.argsort(rho_tmp))
    rho_rank.append(id_rho)

    plt.plot(id_deg, id_rho, '.', color='gray', markersize=10, alpha=0.8)

deg_arr = deg_rank[0]
rho_arr = rho_rank[0]
for i in np.arange(1, len(deg_rank)):
    deg_arr = np.hstack((deg_arr,deg_rank[i]))
    rho_arr = np.hstack((rho_arr,rho_rank[i]))


# 直线拟合与绘制
a1, b1 = optimize.curve_fit(f_1, deg_arr, rho_arr)[0]
x1 = np.arange(0, np.max(deg_arr), 0.1)
y1 = a1 * x1 + b1
plt.plot(x1, y1, "blue", linewidth=2)
#plt.title(r"$c_5$="+str(c5)+", positive correlation", fontsize=15)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r"Ranked $d_{i_s}$", fontsize=20)
plt.ylabel(r"Ranked $r_N$", fontsize=20)
plt.savefig("Figures_sigma_001_1/ranked_global_synch_vs_deg_2.pdf")
print(sst.pearsonr(deg_arr,rho_arr))

C_rank = []
plt.figure(figsize=[6.4 * 1.2, 4.8 * 1.4])
for i in range(n_trial):
    C_list = np.load("data_results_sigma_001_1/WCOs_chimeta_id_" + str(c5) + "_aal2_N94_"+str(i)+".npy")
    C_tmp = C_list[chimera_id[i]]
    deg_tmp = deg_mat[chimera_id[i]]
    id_deg = np.argsort(np.argsort(deg_tmp))
    id_C = np.argsort(np.argsort(C_tmp))
    C_rank.append(id_C)
    plt.plot(id_deg, id_C, '.', color='gray', markersize=10, alpha=0.8)

C_arr = C_rank[0]
for i in np.arange(1, len(C_rank)):
    C_arr = np.hstack((C_arr,C_rank[i]))

# 直线拟合与绘制
a2, b2 = optimize.curve_fit(f_1, deg_arr, C_arr)[0]
x2 = np.arange(0, np.max(deg_arr), 0.1)
y2 = a2 * x2 + b2
#plt.plot(x2, y2, "blue", linewidth=2)
#plt.title(r"$c_5$="+str(c5)+", weak negative correlation", fontsize=15)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r"Ranked $d_{i_s}$", fontsize=20)
plt.ylabel(r"Ranked $\Gamma$", fontsize=20)
plt.savefig("Figures_sigma_001_1/ranked_chimera_id_vs_deg_2.pdf")
print(sst.pearsonr(deg_arr,C_arr))


"""
c5 = 330
L_idxy = np.load("data_results_sigma_001_1/WCOs_index_lower_threshold_c5_" + str(c5)+".npy")
meta_id = np.array([np.arange(94)]*n_trial)
meta_id[L_idxy[0], L_idxy[1]] = -N
meta_id = meta_id < 0

lambda_rank = []
deg_rank = []
plt.figure(figsize=[6.4 * 1.2, 4.8 * 1.4])
for i in range(n_trial):
    lambda_list = np.load("data_results_sigma_001_1/WCOs_chimera_id_"+str(c5)+"_aal2_N94_"+str(i)+".npy")
    lambda_tmp = lambda_list[meta_id[i]]
    deg_tmp = deg_mat[meta_id[i]]
    id_deg = np.argsort(np.argsort(deg_tmp))
    id_lambda = np.argsort(np.argsort(lambda_tmp))
    lambda_rank.append(id_lambda)
    deg_rank.append(id_deg)
    plt.plot(id_deg, id_lambda, '.', color='gray', markersize=10, alpha=0.8)

deg_arr = deg_rank[0]
lambda_arr = lambda_rank[0]
for i in np.arange(1, len(C_rank)):
    deg_arr = np.hstack((deg_arr,deg_rank[i]))
    lambda_arr = np.hstack((lambda_arr,lambda_rank[i]))

# 直线拟合与绘制
a3, b3 = optimize.curve_fit(f_1, deg_arr, lambda_arr)[0]
x3 = np.arange(0, np.max(deg_arr), 0.1)
y3 = a3 * x3 + b3
#plt.plot(x3, y3, "blue", linewidth=2)
plt.title(r"$c_5$="+str(c5)+", weak negative correlation", fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r"Ranked weighted degree", fontsize=15)
plt.ylabel(r"Ranked metasability index", fontsize=15)
#plt.savefig("Figures_sigma_001_1/ranked_metastability_id_vs_deg_1.pdf")
print(sst.pearsonr(deg_arr,lambda_arr))

plt.show()


time_terminal = time.time()
print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")
"""