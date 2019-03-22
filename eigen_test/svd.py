import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import csv
from ss_mod import m
from pyomo.environ import *
from jacobi_reader import get_jacobian

# Constructing constraint gradient 

filename = 'jacobi_debug.in'

jacobian,idx2state,idx2eqn = get_jacobian(filename)
nz = len(jacobian)

lam, evec = la.eig(jacobian)
N,s,LT = la.svd(jacobian)
print(s)
#print(LT)
# i^th column of evec is the evector of i^th eigenvalue 

lam_index = np.arange(1,lam.size+1,1)
lam_real = np.real(lam)
np.savetxt('evals_real.txt',np.transpose(np.stack([lam_index,lam_real])), fmt = '%1.2f')
np.savetxt('evals.txt',np.transpose(np.stack([lam_index,lam])),fmt = '%1.2f')

np.savetxt('jacobian.txt',jacobian)

pos_eval = []
for i in range(0,lam.size):
    if np.real(lam[i])>0:
        pos_eval.append(i)

# Express states in basis of right singular vectors
x_1  = np.zeros(nz)
x_2  = np.zeros(nz)
x_3  = np.zeros(nz)
x_4  = np.zeros(nz)
x_5  = np.zeros(nz)
x_6  = np.zeros(nz)
x_7  = np.zeros(nz)
x_8  = np.zeros(nz)
x_9  = np.zeros(nz)
x_10 = np.zeros(nz)
x_11 = np.zeros(nz)

x_1[0]   = 1
x_2[1]   = 1
x_3[2]   = 1
x_4[3]   = 1
x_5[4]   = 1
x_6[5]   = 1
x_7[6]   = 1
x_8[7]   = 1
x_9[8]   = 1
x_10[9]  = 1
x_11[10] = 1

x_1l  = LT.dot(x_1)
x_2l  = LT.dot(x_2)
x_3l  = LT.dot(x_3)
x_4l  = LT.dot(x_4)
x_5l  = LT.dot(x_5)
x_6l  = LT.dot(x_6)
x_7l  = LT.dot(x_7)
x_8l  = LT.dot(x_8)
x_9l  = LT.dot(x_9)
x_10l = LT.dot(x_10)
x_11l = LT.dot(x_11)


with open('jac_svd.txt','w') as f:
    f.write('sing. vals:')
    f.write('\t')
    for i in range(1,nz+1):
        f.write('%1.3e'%s[i-1]+'\t')
    f.write('\n\n')
    f.write('state')
    f.write('\t\t')
    f.write('contribution')
    f.write('\n')
    for i in range(1,nz+1):
        f.write(idx2state[i])
        x = eval('x_{0}l'.format(str(i)))
        print(x)
        for j in range(0,nz):
            f.write('\t')
            f.write('%1.3e'%(x[j]))
        f.write('\n')



#print(pos_eval)

#def idx2state(idx):
#    # idx is the index, according to numpy, in an eigenvector (array)
#    M1 = 41
#    M2 = M1+41
#    x1 = M2+82
#    x2 = x1+82
#    if idx >= 0 and idx < M1:
#        i = idx+1
#        return 'M1[%i]'%i
#    if idx >= M1 and idx < M2:
#        return 'M2[%i]'%(idx-M1+1)
#    if idx >= M2 and idx < x1:
#        return 'x1[%i]'%(idx-M2+1)
#    if idx >= x1 and idx < nz:
#        return 'x2[%i]'%(idx-x1+1)
#    if idx >= nz or idx < 0:
#        raise ValueError('idx2state function must receive a valid index for a differential variable')

'''

idx2state=  {}
nM1 = 41
nM2 = 41
nx1 = 82
nx2 = 82
idx = 0
for coord in range(0,nM1):
    idx2state[coord+idx] = 'M1[%i]'%(coord+1)
idx = idx + nM1
for coord in range(0,nM2):
    idx2state[coord+idx] = 'M2[%i]'%(coord+1)
idx = idx + nM2
for coord in range(0,nx1):
    if coord & 1:
        idx2state[coord+idx] = 'x1[%i,2]'%((coord+1)/2)
    else:
        idx2state[coord+idx] = 'x1[%i,1]'%(coord/2+1)
idx = idx + nx1
for coord in range(0,nx2):
    if coord & 1:
        idx2state[coord+idx] = 'x2[%i,2]'%((coord+1)/2)
    else:
        idx2state[coord+idx] = 'x2[%i,1]'%(coord/2+1)

#print(idx2state)

# Really should have idx2state dictionary, so it can be given to ol_run file to interpret eigenvector
# ... could just pass over states_impacted[] ...
# No, because I need to assign states to vector coordinates
        
states_impacted = {}
for ev_num in pos_eval:
    states_impacted[ev_num] = []
    ev = evec[:,ev_num]
    r_ev = np.real(ev)
    for i in range(0,lam.size):
        if r_ev[i] > 1e-7 or r_ev[i] < -1e-7:
            states_impacted[ev_num].append(idx2state[i])
 
#print(states_impacted[2])
#print(evec[:,2])

sum_ptb = np.zeros(nz)
for j in range(0,nz):
    sum_ptb[j] = sum(evec[i,j] \
        for i in range(int(var_locator_z['M1'][0])-2,int(var_locator_z['M1'][1])))

#sum_ptb = sum(evec[i,2] \
#        for i in range(int(var_locator_z['M1'][0]-1),int(var_locator_z['M1'][1])))
#print(sum_ptb)



#lam1,evec1 = la.eig(dfdz) 
#lam1_index = np.arange(1,lam1.size+1,1)
#lam1_real = np.real(lam1)
#fig2,ax2 = plt.subplots()
#eval1_spec = ax2.bar(lam1_index,lam1_real)
#plt.show()

#print(lam)
#print(evec)

evec_inv_T = np.transpose(la.inv(evec))
UPSR = evec*evec_inv_T
#np.savetxt('UPSR.txt',np.real(UPSR))

#print(UPSR)
'''
fig1,ax1 = plt.subplots()
eval_spec = ax1.bar(lam_index,lam_real)
plt.ylim(-70,70)
fig1.suptitle('Eigenvalues of Jacobian')
fig1.savefig('eigenvalues.png')

lam_sort = np.sort(np.abs(lam_real))
'''
fig2,ax2 = plt.subplots()
ax2.set_yscale('log')
eval_spec_sort = ax2.bar(lam_index,lam_sort)
fig2.suptitle('Abs of Eigenvalues, sorted')
fig2.savefig('eigenvalues_sort.png')

def state2UPSR(var,no):
# no refers to the index of the variable.
# for example, to access the row corresponding 
# to feed tray holdup, call ('M1',21)
    start = var_locator_z[var][0]
    end = var_locator_z[var][1]
    if var in Mvar:
        start = start - 1
        end = end + 1
        # ^ account for fact that reboiler and condenser holdups don't show up
    if no < 1 or no > 1 + end-start:
        raise ValueError('No. must be a valid index')
    return UPSR[ int(start + (no-1) - 1), :] 

'''
