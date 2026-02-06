from arc import *                 # 导入 ARC（Alkali Rydberg Calculator）
import math 
from math import *
import numpy as np
from scipy.io import savemat
from dataclasses import dataclass
import matplotlib.pyplot as plt
import mpmath

spectrum_point = np.zeros([9157,2])
spectrum_value = np.linspace(450, 650, 9157)
spectrum_point[:,0] = spectrum_value
mpmath.mp.dps = 100
pi = math.pi
atom = Cesium()
c = 299792458
n_min = 4
n_max = 50
level = np.zeros([10000,4])
NUM = 0
P=[1e-6,2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,9e-6,1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5,1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,7e-2,8e-2,9e-2,1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,1e0,2e0,3e0,4e0,5e0,6e0,7e0,8e0,9e0]
# P=[1000]
for n2 in range(n_min, n_max):
    for l2 in range(0, n2):
        if l2 != 0:
            for j2 in [l2-0.5, l2+0.5]:
                transition_wavelength = atom.getTransitionWavelength(6, 0, 0.5, n2, l2, j2) * 1e9
                # print("n2=%d, l2=%d, j2=%.1f" % (n2, l2, j2), end=' ')
                # print("transition_wavelength = %.10f nm" % transition_wavelength, end=' \n')
                level[NUM,0] = n2
                level[NUM,1] = l2
                level[NUM,2] = j2
                level[NUM,3] = transition_wavelength
                NUM = NUM + 1
        else:
            transition_wavelength = atom.getTransitionWavelength(6, 0, 0.5, n2, 0, 0.5) * 1e9
            # print("n2=%d, l2=%d, j2=%.1f" % (n2, 0, 0.5), end=' ')
            # print("transition_wavelength = %.10f nm" % transition_wavelength, end=' \n')
            level[NUM,0] = n2
            level[NUM,1] = 0
            level[NUM,2] = 0.5
            level[NUM,3] = transition_wavelength
            NUM = NUM + 1
Cs_level = level[(-level[:, 3]).argsort()]
Cs_level = np.delete(Cs_level, np.arange(NUM-5, NUM), axis=0)
Cs_level = np.delete(Cs_level, 12, axis=0)
Cs_level = np.delete(Cs_level, 6, axis=0)
Cs_level = np.delete(Cs_level, 1, axis=0)

NUM = 2300
Frequence_level = np.zeros(NUM)
for index in range(0, NUM):
    Frequence_level[index] = atom.getTransitionFrequency(6, 0, 0.5, int(Cs_level[index,0]), int(Cs_level[index,1]), Cs_level[index,2])

dipole_element = np.zeros([NUM,NUM])
for index1 in range(0, NUM):
    for index2 in range(0, NUM):
        if (np.abs(int(Cs_level[index1,1]) - int(Cs_level[index2,1]))>1) | (np.abs(Cs_level[index1,2] - Cs_level[index2,2])>1):
            dipole_element[index1, index2] = 0
        elif Cs_level[index1,2] < Cs_level[index2,2]:
            dipole_element[index1, index2] = np.abs(atom.getReducedMatrixElementJ(int(Cs_level[index1,0]), int(Cs_level[index1,1]), floor(Cs_level[index1,2])+0.5, int(Cs_level[index2,0]), int(Cs_level[index2,1]), floor(Cs_level[index2,2])+0.5))
        else:
            dipole_element[index1, index2] = np.abs(atom.getReducedMatrixElementJ(int(Cs_level[index2,0]), int(Cs_level[index2,1]), floor(Cs_level[index2,2])+0.5, int(Cs_level[index1,0]), int(Cs_level[index1,1]), floor(Cs_level[index1,2])+0.5))
        # print("num=%d -> %d" % (index1, index2), end=' ')

lifetime = np.zeros(NUM)
for index in range(0, NUM):
    if Cs_level[index,1] < 2:
        lifetime[index] = atom.getStateLifetime(int(Cs_level[index,0]), int(Cs_level[index,1]), Cs_level[index,2]) * 1
    else:
        lifetime[index] = atom.getStateLifetime(int(Cs_level[index,0]), int(Cs_level[index,1]), Cs_level[index,2]) * 1
        
A = np.zeros([NUM,NUM])
for index1 in range(0, NUM):
    for index2 in range(0, index1):
        # A[index1,index2] = 1 / lifetime[index1] * dipole_element[index1, index2]**2 * (Frequence_level[index1] - Frequence_level[index2])**3 / np.sum((-1*Frequence_level[0:index1] + Frequence_level[index1])**3 * dipole_element[index1, 0:index1]**2)
        A[index1,index2] = atom.getTransitionRate(int(Cs_level[index1,0]), int(Cs_level[index1,1]), Cs_level[index1,2], int(Cs_level[index2,0]), int(Cs_level[index2,1]), Cs_level[index2,2], temperature=300.0)

# ##############################################################
###############################################################
##6S1/2 -> 6P3/2 -> 24D3/2 -> 26P1/2
step1 = 2
step2 = 728
step3 = 838
wavelength1 = 497.925
wavelength2 = 595.936
spectrum_point[:,1] = 0
d1 = dipole_element[0,step1]
d2 = dipole_element[step1,step2]
d3 = dipole_element[step2,step3]

p1 = 10   #mw
p2 = 10  #mw
p3 = 0e-2 #mw
W1 = 1e6
W2 = (d2/d1)**2*p2/p1*W1
W3 = (d3/d1)**2*p3/p1*W1

M = np.zeros([step3+1,step3+1])
for index1 in range(0, step3):
    for index2 in range(0, step3+1):
        if index2 < index1:
            M[index1, index2] = -A[index2,index1]
        elif index2 == index1:
            M[index1, index2] = -1 * np.sum(A[index1,0:index2])
        else:
            M[index1, index2] = A[index2,index1]
M[step3,:] = 1
M[0,0] = M[0,0] - W1
M[0,step1] = M[0,step1] + W1
M[step1,0] = M[step1,0] + W1
M[step1,step1] = M[step1,step1] - W2 - W1
M[step1,step2] = M[step1,step2] + W2
M[step2,step1] = M[step2,step1] + W2
M[step2,step2] = M[step2,step2] - W3 - W2
M[step2,step3] = M[step2,step3] + W3
M_inv = np.linalg.inv(M)
N = 6.02e23
N_BJ = M_inv[:,step3]/M_inv[0,step3] * N

#M[98,0] = M[98,0] + W_1

sp_wavelength = np.zeros([step3+1,step3+1])
for index1 in range(0, step3+1):
    for index2 in range(0, index1):
        sp_wavelength[index1, index2] = c / (Frequence_level[index1] - Frequence_level[index2] + 1e-6) * 1e9

sp_I = np.zeros([step3+1,step3+1])
for index1 in range(0, step3+1):
    for index2 in range(0, index1):
        sp_I[index1, index2] = (N_BJ[index1]) * A[index1, index2] * (Frequence_level[index1] - Frequence_level[index2]) * 6.65e-34

SP = np.zeros([10000,8])
num = 0
for index1 in range(0, step3+1):
    for index2 in range(0, index1):
        if (sp_I[index1, index2] > 0) & (sp_wavelength[index1, index2] < 650) & (sp_wavelength[index1, index2] > 450):
            # print("index1 = %.5f, index2 = %.5f, wavelength = %.5f" % (index1, index2, sp_wavelength[index1, index2]), end='\n ')
            SP[num,0] = sp_I[index1, index2]
            SP[num,1] = sp_wavelength[index1, index2]
            SP[num,2] = Cs_level[index1,0]
            SP[num,3] = Cs_level[index1,1]
            SP[num,4] = Cs_level[index1,2]
            SP[num,5] = Cs_level[index2,0]
            SP[num,6] = Cs_level[index2,1]
            SP[num,7] = Cs_level[index2,2]
            num = num + 1
SP_off = SP[:,:].copy()   
 
spectrum_point[:,1] = 0 
for i in range(0,num):
    for j in range(0,9156):
        if (SP_off[i,1] <= spectrum_point[j+1,0]) & (SP_off[i,1] >= spectrum_point[j,0]):
            spectrum_point[j,1] = SP_off[i,0] + spectrum_point[j,1]
plt.vlines(spectrum_point[:,0], ymin=0, ymax=spectrum_point[:,1], colors='red', linewidth=1)
# 添加标题和坐标轴标签
plt.title("THzoff")  # 图表标题
plt.xlabel("wavelength (nm)")      # X轴标签
plt.ylabel("Idensity")      # Y轴标签
# 显示图形
plt.show()      
    #-----------------------------------#
D = np.zeros([63])
S = np.zeros([63])
AA = np.zeros([63,3])
D_1 = np.zeros([63])
S_1 = np.zeros([63])
BB = np.zeros([63,3])
index_p = 0
for p3 in P:    
    W3 = (d3/d1)**2*p3/p1*W1    
    M = np.zeros([step3+1,step3+1])
    for index1 in range(0, step3):
        for index2 in range(0, step3+1):
            if index2 < index1:
                M[index1, index2] = 0
            elif index2 == index1:
                M[index1, index2] = -1 * np.sum(A[index1,0:index2])
            else:
                M[index1, index2] = A[index2,index1]
    M[step3,:] = 1
    M[0,0] = M[0,0] - W1
    M[0,step1] = M[0,step1] + W1
    M[step1,0] = M[step1,0] + W1
    M[step1,step1] = M[step1,step1] - W2 - W1
    M[step1,step2] = M[step1,step2] + W2
    M[step2,step1] = M[step2,step1] + W2
    M[step2,step2] = M[step2,step2] - W3 - W2
    M[step2,step3] = M[step2,step3] + W3
    M_inv = np.linalg.inv(M)
    N = 6.02e23
    N_BJ = M_inv[:,step3]/M_inv[0,step3] * N
    
    #M[98,0] = M[98,0] + W_1
    
    sp_wavelength = np.zeros([step3+1,step3+1])
    for index1 in range(0, step3+1):
        for index2 in range(0, index1):
            sp_wavelength[index1, index2] = c / (Frequence_level[index1] - Frequence_level[index2] + 1e-6) * 1e9
    
    sp_I = np.zeros([step3+1,step3+1])
    for index1 in range(0, step3+1):
        for index2 in range(0, index1):
            sp_I[index1, index2] = (N_BJ[index1]) * A[index1, index2] * (Frequence_level[index1] - Frequence_level[index2]) * 6.65e-34
    
    SP = np.zeros([10000,8])
    num = 0
    for index1 in range(0, step3+1):
        for index2 in range(0, index1):
            if (sp_I[index1, index2] > 0) & (sp_wavelength[index1, index2] < 650) & (sp_wavelength[index1, index2] > 450):
                # print("index1 = %.5f, index2 = %.5f, wavelength = %.5f" % (index1, index2, sp_wavelength[index1, index2]), end='\n ')
                SP[num,0] = sp_I[index1, index2]
                SP[num,1] = sp_wavelength[index1, index2]
                SP[num,2] = Cs_level[index1,0]
                SP[num,3] = Cs_level[index1,1]
                SP[num,4] = Cs_level[index1,2]
                SP[num,5] = Cs_level[index2,0]
                SP[num,6] = Cs_level[index2,1]
                SP[num,7] = Cs_level[index2,2]
                num = num + 1
    SP_on = SP[:,:].copy()    
    
    spectrum_point[:,1] = 0 
    for i in range(0,num):
        for j in range(0,9156):
            if (SP_on[i,1] <= spectrum_point[j+1,0]) & (SP_on[i,1] >= spectrum_point[j,0]):
                spectrum_point[j,1] = SP_on[i,0] + spectrum_point[j,1]
    plt.vlines(spectrum_point[:,0], ymin=0, ymax=spectrum_point[:,1], colors='red', linewidth=1)
    # 添加标题和坐标轴标签
    plt.title("THzon")  # 图表标题
    plt.xlabel("wavelength (nm)")      # X轴标签
    plt.ylabel("Idensity")      # Y轴标签
    # 显示图形
    plt.show()    
    
    Diff = SP_on.copy()
    for i in range(0,num):
        for j in range (i,num):
            if Diff[i,1]==Diff[i+j,1]:
                Diff[i+j,0] = Diff[i+j,0] + Diff[i,0]
                Diff[i,1] = 0
                Diff[i,0] = 0
        
    for i in range(0,num):
        for j in range(0,10000):
            if (Diff[i,1] == SP_off[j,1]):
                Diff[i,0] = Diff[i,0] - SP_off[j,0]
                
    spectrum_point[:,1] = 0 
    for i in range(0,num):
        for j in range(0,9156):
            if (Diff[i,1] <= spectrum_point[j+1,0]) & (Diff[i,1] >= spectrum_point[j,0]):
                spectrum_point[j,1] = Diff[i,0] + spectrum_point[j,1]
    
    # plt.vlines(Diff[:,1], ymin=0, ymax=Diff[:,0], colors='red', linewidth=1)
    # # 添加标题和坐标轴标签
    # plt.title("THzon - THzoff（theory）")  # 图表标题
    # plt.xlabel("wavelength (nm)")      # X轴标签
    # plt.ylabel("Idensity")      # Y轴标签
    # # 显示图形
    # plt.show() 
                
    # plt.vlines(spectrum_point[:,0], ymin=0, ymax=spectrum_point[:,1], colors='red', linewidth=1)
    # # 添加标题和坐标轴标签
    # plt.title("THzon - THzoff（experiment）")  # 图表标题
    # plt.xlabel("wavelength (nm)")      # X轴标签
    # plt.ylabel("Idensity")      # Y轴标签
    # # 显示图形
    # plt.show()   
    
    for i in range(0,num):
        if (np.abs(Diff[i,1] - wavelength1)<1e-3):
            D[index_p] = Diff[i,0]
    for i in range(0,num):
        if (np.abs(SP_on[i,1] - wavelength1)<1e-3):
            S[index_p] = SP_on[i,0]        
    for i in range(0,num):
        if (np.abs(Diff[i,1] - wavelength2)<1e-2):
            D_1[index_p] = Diff[i,0]
    for i in range(0,num):
        if (np.abs(SP_on[i,1] - wavelength2)<1e-2):
            S_1[index_p] = SP_on[i,0]        
    index_p = index_p + 1
    print("pthz=%f" % (p3), end=' ')
# #########################################################################
# plt.plot(np.log10(P[:]), np.log10(S[:]))
# plt.show()   
# plt.plot(np.log10(P[:])-np.log10(P[0]), np.log10(np.abs(D[:]))-np.log10(np.abs(D[0])))
# plt.show()   
AA[:,0] = P
AA[:,1] = S
AA[:,2] = D
BB[:,0] = P
BB[:,1] = S_1
BB[:,2] = D_1
savemat('30D32_33P12_wavelength1.mat', {'wavelength1_497_925_30D32_33P12': AA})
savemat('30D32_33P12_wavelength2.mat', {'wavelength2_595_936_30D32_33P12': BB})
# savemat('24D32_26P12_thzoff.mat', {'SPthzoff_24D32_26P12': SP_off})
# savemat('24D32_26P12_thzon.mat', {'SPthzon_24D32_26P12': SP_on})