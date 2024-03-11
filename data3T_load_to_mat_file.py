# -*- coding: utf-8 -*-
from typing import TextIO

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.io import savemat
import scipy.stats as sst
import seaborn as sns
import time


time_start = time.time()

#data_list = ['100610',]

f = open('E:/dsi_studio/data_3T/list.txt')
data_id = f.read().splitlines()
f.close()



N = 120
A_count = np.zeros((N, N))
count_max = []
compare_name = []
compare_atlas = []
for i in data_id:
    data_count = sio.loadmat('E:/dsi_studio/data_3T/'+i+'/T1w/Diffusion/whole_brain_AAL2_count.mat')
    A_count += data_count["connectivity"]
    count_max.append(np.max(data_count["connectivity"]))
    if i==data_id[0]:
        name_0 = data_count['name']
        atlas_0 = data_count['atlas']
    else:
        compare_name.append(np.allclose(name_0, data_count['name']))
        compare_atlas.append(np.allclose(atlas_0, data_count['atlas']))

A_count /= len(data_id)
A_count[A_count<np.max(A_count)*0.001] = 0
data_count['connectivity'] = A_count
savemat('data3T_aal2_connectivity_count.mat', mdict=data_count,  format='4')

N = 94
A_count = A_count[:N, :N]
A_count = A_count/np.sum(A_count/2)
np.save("data3T_aal2_A_count_N94_1.npy", A_count)

voxel_size = 2*10**-3
dist_mat = np.load("dist_mat_aal2.npy")
dist_mat *= voxel_size
dist_mat = dist_mat[:94, :94]

"""
print('compare_name :', compare_name)
print(list(map(int, compare_name)))
print('compare_atlas :', compare_atlas)
print(list(map(int, compare_atlas)))

A_arr = A_count.flatten('C')
dist_arr = dist_mat.flatten('C')
id_A = np.where(A_arr>0)[0]
print(sst.pearsonr(A_arr[id_A],dist_arr[id_A]))

"""
time_terminal = time.time()
print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")


plt.figure()
plt.imshow(A_count, cmap='Blues')
cbar = plt.colorbar()
tick_font_size = 15
cbar.ax.tick_params(labelsize=tick_font_size)
plt.title(r"Connectivity Matrix", fontsize=15)
plt.gca().invert_yaxis()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.savefig("Figures_sigma_001_1/connectivity_matrix_aal2_N94_1.pdf")

plt.figure()
plt.imshow(dist_mat, cmap='Reds')
cbar = plt.colorbar()
tick_font_size = 15
#cbar.set_label(r'spatial distance (m)')
cbar.set_label(r'$d_{ij}$ (m)')
cbar.ax.tick_params(labelsize=tick_font_size)
cbar.ax.figure.axes[-1].yaxis.label.set_size(15)

plt.title(r"Distance Matrix", fontsize=15)
plt.gca().invert_yaxis()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.savefig("Figures_sigma_001_1/distance_matrix_aal2_N94_1.pdf")



plt.show()