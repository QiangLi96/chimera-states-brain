# -*- coding: utf-8 -*-
import numpy as np
import xlrd
from collections import defaultdict
import matplotlib.pyplot as plt
import math
import seaborn as sns
import time

time_start = time.time()


c5 = 330
n_trial = 10

# brain_regions_cognitive_systems_N94.xlsx存储着脑区与认知系统的对应信息
book_1=xlrd.open_workbook('D:/article/time_9_20220118/brain_regions_cognitive_systems_N94.xlsx')
sheet_1=book_1.sheet_by_index(0)
id = sheet_1.col_values(0)  # brain regions index
cs = sheet_1.col_values(2)  # cognitive system names
cs_br = defaultdict(list)  # cognitive systems with brain regions index

for i in range(sheet_1.nrows):
    id = int(sheet_1.cell(i,0).value)-1
    cs_br[str(sheet_1.cell(i,2).value)].append(id)

#cs_names = list(cs_br.keys())
cs_names = ['Som', 'DMN', 'Con', 'VA', 'Lim', 'Oth', 'Vis', 'DA']


s_cog = list(cs_br.values())

cog_tmp = s_cog[3]
s_cog[3] = s_cog[5]
s_cog[5] = cog_tmp

name_tmp = cs_names[3]
cs_names[3] = cs_names[5]
cs_names[5] = name_tmp

cog_tmp = s_cog[4]
s_cog[4] = s_cog[-1]
s_cog[-1] = cog_tmp

name_tmp = cs_names[4]
cs_names[4] = cs_names[-1]
cs_names[-1] = name_tmp

rho_list=[]
C_list=[]
for k in range(len(s_cog)):
    rho_list.append([])
    C_list.append([])

chimera_id = np.load("data_results_sigma_001_1/WCOs_index_chimera_c5_" + str(c5)+".npy")

for i in range(n_trial):
    results_rho = np.load("data_results_sigma_001_1/WCOs_global_synch_" + str(c5) + "_aal2_N94_"+str(i)+".npy")
    C_arr = np.load("data_results_sigma_001_1/WCOs_chimeta_id_" + str(c5) + "_aal2_N94_" + str(i) + ".npy")
    for j in range(len(s_cog)):
        id_tmp = chimera_id[i][s_cog[j]]
        sc_tmp = np.array(s_cog[j])
        sc_tmp = sc_tmp[id_tmp]
        rho_list[j].append(results_rho[sc_tmp])
        C_list[j].append(C_arr[sc_tmp])

rho_mean = []
rho_std = []
C_mean = []
C_std = []

x0 = 0
y0 = 0

for k in range(len(rho_list)):
    rho_tmp = []
    C_tmp = []
    for l in range(len(rho_list[k])):
        rho_tmp += rho_list[k][l].tolist()
        C_tmp += C_list[k][l].tolist()
        x0 += len(rho_list[k][l])
        y0 += len(C_list[k][l])
    rho_mean.append(np.mean(rho_tmp))
    rho_std.append(np.std(rho_tmp))
    C_mean.append(np.mean(C_tmp))
    C_std.append(np.std(C_tmp))

rho_sort = np.argsort(rho_mean)
C_sort = np.argsort(C_mean)
print("global synchrony sort: ", np.array(cs_names)[rho_sort])
print("chimera index sort: ", np.array(cs_names)[C_sort])



time_terminal = time.time()
print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")


plt.figure(figsize=[6.4, 4.8*1.4])
plt.bar(np.arange(len(cs_names)), rho_mean,
        alpha=1.0, width = 0.6, edgecolor = 'white', lw=1)
plt.errorbar(np.arange(len(cs_names)), rho_mean,
             alpha=1.0,yerr=rho_std,fmt='o',ecolor='k',color='k',elinewidth=2,capsize=4)
plt.plot(np.arange(len(cs_names)), rho_mean, color='gray', alpha=0.6)
plt.xticks(np.arange(len(cs_names)), labels=cs_names, fontsize=20, rotation=45)
plt.yticks(fontsize=20)
#plt.xlim(-1, len(cs_names)+1)
#plt.ylim(0, np.max(np.array(rho_mean)+np.array(rho_std))+0.2)
#plt.legend(loc="upper right", fontsize=15)
plt.savefig("Figures_sigma_001_1/global_synch_vs_cs_hist_1.pdf")


plt.figure(figsize=[6.4, 4.8*1.4])
plt.bar(np.arange(len(cs_names)), C_mean,
        alpha=1.0, width = 0.6, edgecolor = 'white', lw=1)
plt.errorbar(np.arange(len(cs_names)),C_mean,
             alpha=1.0,yerr=C_std,fmt='o',ecolor='k',color='k',elinewidth=2,capsize=4)
plt.plot(np.arange(len(cs_names)),C_mean, color='gray', alpha=0.6)
plt.xticks(np.arange(len(cs_names)), labels=cs_names, fontsize=20, rotation=45)
plt.yticks(fontsize=20)
#plt.xlim(-1, len(cs_names)+1)
#plt.ylim(0, np.max(np.array(rho_mean)+np.array(rho_std))+0.2)
#plt.legend(loc="upper right", fontsize=15)
plt.savefig("Figures_sigma_001_1/chimera_id_vs_cs_hist_1.pdf")



"""
plt.figure(figsize=[6.4 * 1.2, 4.8 * 1.2])
plt.bar(np.arange(len(cs_names)), rho_mean,
        alpha=1.0, width = 0.6, facecolor = 'red', edgecolor = 'white', label=r'$r_N$', lw=1)
plt.errorbar(np.arange(len(cs_names)), rho_mean,
             alpha=1.0,yerr=rho_std,fmt='o',ecolor='k',color='k',elinewidth=2,capsize=4)
plt.plot(np.arange(len(cs_names)), rho_mean, color='gray', alpha=0.6)

plt.bar(np.arange(len(cs_names)), C_mean,
        alpha=1.0, width = 0.6, facecolor = 'yellow', edgecolor = 'white', label='$\Gamma$', lw=1)
plt.errorbar(np.arange(len(cs_names)),C_mean,
             alpha=1.0,yerr=C_std,fmt='o',ecolor='k',color='k',elinewidth=2,capsize=4)
plt.plot(np.arange(len(cs_names)),C_mean, color='gray', alpha=0.6)
#plt.title(r"$c_5$="+str(c5), fontsize=15)
plt.xticks(np.arange(len(cs_names)), labels=cs_names, fontsize=15, rotation=45)
plt.yticks(fontsize=15)
#plt.xlim(-1, len(cs_names)+1)
#plt.ylim(0, np.max(np.array(rho_mean)+np.array(rho_std))+0.2)
plt.legend(loc="upper right", fontsize=15)
plt.savefig("Figures_sigma_001_1/global_synch_chimera_id_vs_cs_hist_2.pdf")
"""




