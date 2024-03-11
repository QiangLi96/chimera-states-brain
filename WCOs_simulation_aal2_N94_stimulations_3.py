# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import math
import time

A_mat = np.load("data3T_aal2_A_count_N94_1.npy")

N = len(A_mat)
voxel_size = 2*10**-3


dist_mat = np.load("dist_mat_aal2.npy")
dist_mat *= voxel_size
dist_mat = dist_mat[:94, :94]

td = 10
tau_d = dist_mat / td  # tau_d^{ij}表示j对i的delay，即第i个节点的delay。此处默认 tau_d^{ij} = tau_d^{ji}

t_init = 0
t_end = 150

# Total_steps = 5000  # Compute 1000 grid points
# dt = float(t_end - t_init) / Total_steps
dt = 0.01 * 10 ** (-1)
#Total_steps = int((t_end - t_init) / dt)

E_init = np.ones(N) * 0.1
I_init = np.ones(N) * 0.1

c1 = 16
c2 = 12
c3 = 15
c4 = 3

aE = 1.3
aI = 2
thetaE = 4
thetaI = 3.7
tau = 16

mu = 0
#sigma = 0.00005*100
sigma = 0.001



def fun_SE(x):
    return 1 / (1 + np.exp(-aE * (x - thetaE))) - 1 / (1 + np.exp(aE * thetaE))


# maximal value
S_Em = 1 - 1 / (1 + np.exp(aE * thetaE))


def fun_SI(x):
    return 1 / (1 + np.exp(-aI * (x - thetaI))) - 1 / (1 + np.exp(aI * thetaI))


# maximal value
S_Im = 1 - 1 / (1 + np.exp(aI * thetaI))


def dot_E(t, E_tmp, Ed_tmp, I_tmp, c5_tmp, P_arr):
    E0 = c1 * E_tmp - c2 * I_tmp + c5_tmp * np.sum(np.multiply(A_mat, Ed_tmp), axis=1) + P_arr
    dotE_tmp = -E_tmp + np.multiply((S_Em - E_tmp), fun_SE(E0))
    return dotE_tmp


def dot_I(t, E_tmp, I_tmp, Id_tmp, c6_tmp):
    I0 = c3 * E_tmp - c4 * I_tmp + c6_tmp * np.sum(np.multiply(A_mat, Id_tmp), axis=1)
    dotI_tmp = -I_tmp + np.multiply((S_Im - I_tmp), fun_SI(I0))
    return dotI_tmp


def dW(delta_t):
    """Sample a random number at each call."""
    return np.random.normal(loc=0.0, scale=np.sqrt(delta_t), size=N)


def Ed_Id(t, Es_tmp, Is_tmp):
    t_del = t - tau_d
    i_td = np.around(t_del / dt)
    i_td = i_td.astype(int)
    Ed_mat = np.zeros((N, N))
    Id_mat = np.zeros((N, N))
    for j in range(np.shape(i_td)[1]):
        i_j = i_td[:, j]
        Ed_mat[:, j] = Es_tmp[j, i_j]  # 此处默认 tau_d^{ij} = tau_d^{ji}, Ed_mat 第j列表示第j个节点的delay(相对于其他N-1个节点)
        # Ed_mat 的第i行表示其他N-1个节点的delay(相对于第i个节点)
        Id_mat[:, j] = Is_tmp[j, i_j]  # Id_mat 第j列表示第j个节点的delay(相对其他N-1个节点)
        # Id_mat 的第i行表示其他N-1个节点的delay(相对于第i个节点)
    return Ed_mat, Id_mat


ts = np.arange(t_init, t_end + dt/2, dt)
Total_steps = ts.size
Es = np.zeros((N, Total_steps))
Is = np.zeros((N, Total_steps))

Es[:, 0] = E_init
Is[:, 0] = I_init

k_init = -int(np.round(1 / dt))

def fun_EM(i_trial, c5, i_sti):
    c6 = c5 / 4
    P_arr = np.zeros(N)
    P_arr[i_sti] = 1.15
    for i in range(1, ts.size):
        np.random.seed(i_trial + c5 + i_sti + i + np.random.randint(1000000))
        if ts[i - 1] > np.max(tau_d):
            t = t_init + (i - 1) * dt
            E = Es[:, i - 1]
            I = Is[:, i - 1]
            Ed_arr, Id_arr = Ed_Id(t, Es[:, 0:i], Is[:, 0:i])
            Es[:, i] = E + (dot_E(t, E, Ed_arr, I, c5,P_arr) * dt + sigma * dW(dt))/tau
            Is[:, i] = I + (dot_I(t, E, I, Id_arr, c6) * dt + sigma * dW(dt))/tau
        else:
            t = t_init + (i - 1) * dt
            E = Es[:, i - 1]
            I = Is[:, i - 1]
            Es[:, i] = E + (dot_E(t, E, np.zeros((N, N)), I, c5, P_arr) * dt + sigma * dW(dt)) / tau
            Is[:, i] = I + (dot_I(t, E, I, np.zeros((N, N)), c6) * dt + sigma * dW(dt)) / tau

    phi_t = np.arctan(Is[:, k_init:] / Es[:, k_init:])
    return phi_t


c5 = [1000]
n_trial = 10
c5_Isti = [(x, y, z) for x in range(n_trial) for y in c5 for z in range(N)]
cores = mp.cpu_count()
if __name__ == '__main__':
    time_start = time.time()
    with mp.Pool(cores) as p:
        results = p.starmap(fun_EM, c5_Isti)
    for i in range(n_trial):
        phi_arr = np.array(results[i*N:(i+1)*N])
        np.save("data_simulation_1/WCOs_aal2_N94_stimulations_c5_" + str(c5[0]) + "_phi_arr_"+str(i)+".npy", phi_arr)
    time_terminal = time.time()
    print('totally cost', str("{:.2f}".format(time_terminal - time_start)) + "s")
