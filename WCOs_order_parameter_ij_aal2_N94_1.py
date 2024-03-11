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
i_trial = 2

Es = np.load("WCOs_simulation_Es_aal2_"+str(int(c5))+"_N94_sigma_001_Pi_0_"+str(i_trial)+".npy")
Is = np.load("WCOs_simulation_Is_aal2_"+str(int(c5))+"_N94_sigma_001_Pi_0_"+str(i_trial)+".npy")

N = len(Es)

t_init = 0
t_end = 150
dt = 0.01 * 10 ** (-1)

c6 = c5 / 4

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


ts = np.arange(t_init, t_end + dt/2, dt)
Total_steps = ts.size

ts_mean = ts[-int(np.round(1/dt)):]
k_init = -int(np.round(1/dt))

rho_ijt = np.zeros((len(s_cog),len(s_cog)))
for i in range(len(s_cog)):
    s_i = s_cog[i]
    for j in range(len(s_cog)):
        s_j = s_cog[j]
        for k in np.arange(k_init,0):
            phi_t = np.arctan(Is[:, k] / Es[:, k])
            rho_ijt[i, j] += np.sqrt((np.sum(np.cos(phi_t[s_i])) + np.sum(np.cos(phi_t[s_j]))) ** 2 + \
                            (np.sum(np.sin(phi_t[s_i])) + np.sum(np.sin(phi_t[s_j]))) ** 2) / (len(s_i) + len(s_j))


rho_ij = rho_ijt/(len(ts_mean))



time_terminal = time.time()
print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")


# figure
plt.figure(figsize=[6.4*1.4,4.8*1.5])
# plot heatmap
#ax = sns.heatmap(rho_ij, cmap="viridis", vmin=0.45, vmax=1.0,cbar_kws={'label': r'$r_{\xi_i,\xi_j}$'})
ax = sns.heatmap(rho_ij, vmin=0.58, vmax=1.0,cbar_kws={'label': r'$r_{\xi_i,\xi_j}$'})
cbar = ax.collections[0].colorbar
# here set the labelsize by 15
cbar.ax.tick_params(labelsize=15)
#plt.title(r"$c_5$="+str(c5)+r", $N$="+str(N)+r", $\delta t$="+str(dt)+"\n cognitive systems", fontsize=20)
ax.figure.axes[-1].yaxis.label.set_size(15)
# xticks
ax.xaxis.tick_top()
#xticks_labels = [r"$s_1$", r"$s_2$", r"$s_3$", r"$s_4$", r"$s_5$", r"$s_6$", r"$s_7$", r"$s_8$"]
#plt.xticks(np.arange(len(xticks_labels)) + .5, labels=xticks_labels, fontsize=15)
plt.xticks(np.arange(len(cs_names)) + .5, labels=cs_names, fontsize=20, rotation=45)
# yticks
#ax.yaxis.tick_left()
#yticks_labels = [r"$s_1$", r"$s_2$", r"$s_3$", r"$s_4$", r"$s_5$", r"$s_6$", r"$s_7$", r"$s_8$"]
#plt.yticks(np.arange(len(yticks_labels)) + .5, labels=yticks_labels, fontsize=15)
plt.yticks(np.arange(len(cs_names)) + .5, labels=cs_names, fontsize=20, rotation=45)
# axis labels
#plt.ylabel(r"cognitive systems", fontsize=15)

plt.savefig("Figures_sigma_001_1/order_parameter_ij_c5_"+str(int(c5))+"_Pi_0_N94_2.pdf")

plt.show()

