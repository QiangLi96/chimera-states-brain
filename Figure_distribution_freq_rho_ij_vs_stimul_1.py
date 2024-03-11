# -*- coding: utf-8 -*-
import numpy as np
import xlrd
from collections import defaultdict
import matplotlib.pyplot as plt
import math
import seaborn as sns
import time

time_start = time.time()



N=94
c5 = 330
n_trial = 10

rho_Th = 0.8

num = 64

rho_ij = np.load("data_results_sigma_001_1/WCOs_order_parameter_ij_stimulation_c5_" + str(c5)+".npy")

Nu_list = []
Nl_list = []

for i in range(n_trial*N):
    rho_arr = rho_ij[i]
    rho_arr = rho_arr.flatten('C')
    Nu_list.append(len(rho_arr[rho_arr>=rho_Th]))
    Nl_list.append(len(rho_arr[rho_arr<rho_Th]))

Nu_tmp = []
Nl_tmp = []
for j in np.arange(0,n_trial*N, N):
    Nu_tmp.append(np.array(Nu_list[j:j+N]))
    Nl_tmp.append(np.array(Nl_list[j:j+N]))

Nu_tmp = np.array(Nu_tmp)
Nl_tmp = np.array(Nl_tmp)

U_idx, U_idy = np.where(Nu_tmp == num)
U_idxy = np.vstack((U_idx,U_idy))

L_idx, L_idy = np.where(Nl_tmp == num)
L_idxy = np.vstack((L_idx, L_idy))

chimera_id = np.array([np.arange(94)]*n_trial)

chimera_id[U_idx,U_idy] = -N
chimera_id[L_idx, L_idy] = -N
chimera_id = chimera_id >= 0


#np.save("data_results_sigma_001_1/WCOs_index_chimera_c5_" + str(c5)+".npy", chimera_id)
#np.save("data_results_sigma_001_1/WCOs_index_upper_threshold_c5_" + str(c5)+".npy", U_idxy)
#np.save("data_results_sigma_001_1/WCOs_index_lower_threshold_c5_" + str(c5)+".npy", L_idxy)



Nu_tmp[U_idx,U_idy] = - num
Nl_tmp[L_idx,L_idy] = -num

Nu_arr = []
Nl_arr = []
for j in range(N):
    up_tmp = Nu_tmp[:,j]
    Nu_arr.append(up_tmp[up_tmp>0].tolist())
    low_tmp = Nl_tmp[:,j]
    Nl_arr.append(low_tmp[low_tmp>0].tolist())

"""
x0 = 0
y0 = 0
list0 = []
for k in range(len(Nu_arr)):
    x0 += len(Nu_arr[k])
    y0 += len(Nl_arr[k])
    list0.append(len(Nu_arr[k])==len(Nl_arr[k]))
"""

up_mean = []
up_min = []
up_max = []

low_mean = []
low_min = []
low_max = []

i_sti = 1
sti_list = []
for k in range(len(Nu_arr)):
    if len(Nu_arr[k])>0:
        up_mean.append(np.mean(Nu_arr[k])/num)
        up_min.append(np.min(Nu_arr[k])/num)
        up_max.append(np.max(Nu_arr[k])/num)

        low_mean.append(np.mean(Nl_arr[k])/num)
        low_min.append(np.min(Nl_arr[k])/num)
        low_max.append(np.max(Nl_arr[k])/num)

        sti_list.append(i_sti)

    i_sti += 1

no_sti = np.arange(1,94+1).tolist()




time_terminal = time.time()
print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")

"""
plt.figure(1)
plt.imshow(Nu_tmp/64, vmin=-1, vmax=1)
plt.colorbar()
plt.figure(2)
plt.imshow(Nl_tmp/64, vmin=-1, vmax=1)
plt.colorbar()
plt.figure(3)
plt.imshow(chimera_id/N)
plt.colorbar()
"""

plt.figure(figsize=[6.4*1.2, 4.8*1.4])

current_palette = sns.color_palette("Set1")
sns.set_palette(current_palette)

plt.plot(sti_list,up_mean,'-o', linewidth=2, label=r'$\rho_{s_i,s_j} \geq \rho_{Th}$', alpha=0.8)
#plt.plot(sti_list,low_mean, '-o', label=r'$\rho_{s_i,s_j} < \rho_{Th}$', alpha=0.8)
plt.fill_between(sti_list, up_min, up_max,color=sns.xkcd_rgb['grapefruit'], alpha=0.2)
#plt.fill_between(sti_list, low_min, low_max,color=sns.xkcd_rgb['sky blue'],alpha=0.2)
#plt.plot(sti_list, np.ones(len(sti_list)), '-.', color='gray', alpha=0.6)
#plt.plot(sti_list+[-2, -1, 0], np.zeros(len(sti_list)+3), '-', color='gray', alpha=0.6)
plt.plot(np.arange(-2,np.max(sti_list),0.1), np.zeros(len(np.arange(-2,np.max(sti_list),0.1))), '--', color='gray', alpha=0.6)
#plt.fill_between(mu_range, H_simulation_average-H_simulation_std, H_simulation_average+H_simulation_std,color=sns.xkcd_rgb['grapefruit'],alpha=0.2)
#plt.fill_between(mu_range, H_simulation_average-H_analytic_std, H_simulation_average+H_analytic_std,color=sns.xkcd_rgb['sky blue'],alpha=0.2)
#plt.legend(fontsize=15,loc='upper right')
plt.xlabel(r'$i_s$',fontsize=20)
plt.ylabel(r'$P$',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(-2)
plt.savefig('Figures_sigma_001_1\distribution_proba_rho_ij_vs_stimul_1.pdf')



plt.figure(figsize=[6.4*1.1, 4.8*1.1])
plt.plot(sti_list,up_mean,'-', linewidth=2, label=r'$\rho_{s_i,s_j} \geq \rho_{Th}$', alpha=0.8)
for k in range(len(Nu_arr)):
    if len(Nu_arr[k])>0:
        plt.plot((k+1)*np.ones(len(Nu_arr[k])),np.array(Nu_arr[k])/num, 'o', color='black', alpha=0.2)
plt.xlabel(r'Stimulated node index',fontsize=15)
plt.ylabel(r'Probability',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)


plt.show()


"""
plt.figure()
plt.plot(sti_list,up_mean, 'o', color='red',label='$H$(simulation)')
plt.plot(sti_list,up_min, 'o', color='blue')
plt.plot(sti_list,up_max, 'o', color='yellow')
plt.plot(sti_list,low_mean, '*', color='red', label='$H$(analytic)')
plt.plot(sti_list,low_min, '*', color='blue')
plt.plot(sti_list,low_max, '*', color='yellow')
#plt.fill_between(sti_list, up_min, up_max,color=sns.xkcd_rgb['grapefruit'],alpha=0.2)
#plt.fill_between(sti_list, low_min, low_max,color=sns.xkcd_rgb['sky blue'],alpha=0.2)
plt.plot(sti_list, np.ones(len(sti_list)), '-.', color='gray', alpha=0.6)
plt.plot(sti_list, np.zeros(len(sti_list)), '--', color='gray', alpha=0.6)
plt.legend(fontsize=15)
plt.xlabel(r'simulated node index',fontsize=15)
plt.ylabel(r'frequency',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.savefig('F:\Latex_work\work_2_2019_5_22\Figures_of_article\Fig__H1_Energy_expended_perturbation_range_omega_theta.pdf')
#plt.savefig('F:\Latex_work\work_2_2019_5_22\Figures_of_article\Fig__H1_Energy_expended_perturbation_range_omega_theta.png')
plt.show()
"""