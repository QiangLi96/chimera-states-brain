# -*- coding: utf-8 -*-
import numpy as np
import xlrd
from collections import defaultdict
import multiprocessing as mp
import matplotlib.pyplot as plt
import math
import seaborn as sns
import time

time_start = time.time()



N=94


t_init = 0
t_end = 150
dt = 0.01 * 10 ** (-1)



# brain_regions_cognitive_systems_N94.xlsx存储着脑区与认知系统的对应信息
book_1=xlrd.open_workbook('D:/article/time_9_20220118/brain_regions_cognitive_systems_N94.xlsx')
sheet_1=book_1.sheet_by_index(0)
id = sheet_1.col_values(0)  # brain regions index
cs = sheet_1.col_values(2)  # cognitive system names
cs_br = defaultdict(list)  # cognitive systems with brain regions index

for i in range(sheet_1.nrows):
    id = int(sheet_1.cell(i,0).value)-1
    cs_br[str(sheet_1.cell(i,2).value)].append(id)

cs_names = ['Som', 'DMN', 'Con', 'VA', 'Lim', 'Oth', 'Vis', 'DA']

s_cog = list(cs_br.values())
cog_tmp = s_cog[3]
s_cog[3] = s_cog[-1]
s_cog[-1] = cog_tmp

name_tmp = cs_names[3]
cs_names[3] = cs_names[-1]
cs_names[-1] = name_tmp

cog_tmp = s_cog[4]
s_cog[4] = s_cog[5]
s_cog[5] = cog_tmp

name_tmp = cs_names[4]
cs_names[4] = cs_names[5]
cs_names[5] = name_tmp

ts = np.arange(t_init, t_end + dt/2, dt)
Total_steps = ts.size

ts_mean = ts[-int(np.round(1/dt)):]
k_init = -int(np.round(1/dt))


def fun_rho(phi_tmp):
    rho_ijt = np.zeros((len(s_cog),len(s_cog)))
    for i in range(len(s_cog)):
        s_i = s_cog[i]
        for j in range(len(s_cog)):
            s_j = s_cog[j]
            for k in np.arange(k_init,0):
                phi_t = phi_tmp[:,k]
                rho_ijt[i, j] += np.sqrt((np.sum(np.cos(phi_t[s_i])) + np.sum(np.cos(phi_t[s_j]))) ** 2 + \
                                (np.sum(np.sin(phi_t[s_i])) + np.sum(np.sin(phi_t[s_j]))) ** 2) / (len(s_i) + len(s_j))
    rho_ij = rho_ijt/(len(ts_mean))
    return rho_ij



c5 = 1000
phi_list = []
n_trial = 10
for i in range(n_trial):
    phi_arr = np.load("data_simulation_sigma_001_1/WCOs_aal2_N94_stimulations_c5_"+str(c5)+"_phi_arr_"+str(i)+".npy")
    for j in range(N):
        phi_list.append(phi_arr[j])

cores = mp.cpu_count()
if __name__ == '__main__':
    time_start = time.time()
    with mp.Pool(cores) as p:
        results = p.map(fun_rho, phi_list)

    np.save("data_results_sigma_001_1/WCOs_order_parameter_ij_stimulation_c5_" + str(c5)+".npy", results)
    time_terminal = time.time()
    print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")




