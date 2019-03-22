import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import csv
from ss_mod import m
from pyomo.environ import *

# Constructing constraint gradient 

def get_jacobian(filename):

    jac_file = open(filename,newline='')
    header = jac_file.readline()
    h1 = ''
    h2 = ''
    i=0
    while header[i] != '\t' and header[i] != ' ' and header[i] != '\n':
        h1 += header[i]
        i = i+1
    i = i+1
    while header[i] != '\t' and header[i] != ' ' and header[i] != '\n':
        h2 += header[i]
        i = i+1
    nrow = int(h1)
    ncol = int(h2)
    
    const_jac = np.zeros((nrow,ncol))
    
    reader = csv.reader(jac_file,delimiter='\t')
    
    # put the values from the jacobian file into 
    # their correct places in the jacobian matrix
    for row in reader:
        const_jac[int(row[0])-1,int(row[1])-1] = float(row[2])
    
    #print(const_jac.shape)
    row_sum = np.sum(const_jac,axis=1)
    col_sum = np.sum(const_jac,axis=0)
    
    np.savetxt('row_sum.txt',row_sum)
    np.savetxt('col_sum.txt',col_sum)
    
    np.savetxt('jacobian.txt',const_jac)
    
    rowfilename = 'ss_mod.row'
    colfilename = 'ss_mod.col'
    
    # Parsing row and column files
    
    mbal = ['const11','const12','const15','const18','const20','const32','const33','const36','const39','const41']
    cbal = ['const13','const14','const16','const19','const21','const34','const35','const37','const40','const42']
    vle  = ['const1','const2','const3','const22','const23','const24']
    vflw = ['const4','const5','const25','const26'] 
    lflw = ['const6','const7','const8','const27','const28','const29']
    
    # Differential constraints
    f_list = ['const_E1','const_E2','const_E3','const_S7','const_S8','const_S9','const_S10',
            'const_B7','const_C6','const_C11']
    # Algebraic constraints
    g_list = []
    
    for const in m.component_objects(Constraint):
        if const.name not in f_list:
            g_list.append(const.name)
    
    # Differential variables
    z_list = ['xE_IPA','xS_CO2','TB_tb_ave','TC_tb_ave','TC_sh_ave']
    
    # Algebraic variables
    y_list = []
    
    for var in m.component_objects(Var):
        if var.name not in z_list:
            y_list.append(var.name)
    
    const_list = f_list+g_list
    var_list   = y_list+z_list
    
    
    
    # Maps variables/constraints to arrays holding their start and end line (row/col in jacobian)
    # ^ this seems kinda janky... is there a better way to do it?
    # Note: assumes that indexed constraints are contiguous...
    eqn_locator_jac = {}
    var_locator_jac = {}
    for eqn in const_list: eqn_locator_jac[eqn] = np.zeros(2) 
    for var in var_list:   var_locator_jac[var] = np.zeros(2)
    
    rowfile = open(rowfilename,newline='')
    colfile = open(colfilename,newline='')
    
    # Function to get name of variable/constraint from line of row/col file
    # this also seems super janky
    def get_name(entry):
        i = 0
        name = ''
        while entry[i] != '[' and entry[i] != '\n' and entry[i] != '': 
            name += entry[i]
            i = i+1
        return name    
    
    # Step through lines of both files, recording start and end lines of each equation/variable
    i = 1
    prev = ''
    line = rowfile.readline()
    while line != '': 
        name = get_name(line) 
        if name != prev:
            if name in const_list:
                eqn_locator_jac[name][0] = int(i)
            if prev in const_list:
                eqn_locator_jac[prev][1] = int(i-1)
        prev = name
        line = rowfile.readline()
        i = i+1
    if prev in const_list:
        eqn_locator_jac[prev][1] = int(i-1)
    
        
    
    i = 1
    prev = ''
    line = colfile.readline()
    while line != '': 
        name = get_name(line) 
        if name != prev:
            # dumb hacky fix to problem of possible non-contiguity of variables in .col file
            # definitely needs to be made more thorough/robust
            if name in var_list:
                if var_locator_jac[name][0]==0:
                    var_locator_jac[name][0] = int(i)
            if prev in var_list:
                if var_locator_jac[prev][1]==0:
                    var_locator_jac[prev][1] = int(i-1)
        prev = name
        line = colfile.readline()
        i = i+1
    if prev in var_list:
        if var_locator_jac[prev][1]==0:
            var_locator_jac[prev][1] = int(i-1)
    
    #print(var_locator_jac)
    #print(eqn_locator_jac)
    
    # Construct matrices df,g/dz,y
    
    
    eqn_locator_f = {}
    eqn_locator_g = {}
    var_locator_z = {}
    var_locator_y = {}

    idx2state = {}
    idx2eqn   = {}
    
    # Build index-finder for differential variables in f
    
    # start i at 1 here for consistency with jacobian locator, 
    # which refers to lines in row and col files
    #i = 1
    #for const in mbal:
    #    j = i + eqn_locator_jac[const][1]-eqn_locator_jac[const][0]
    #    eqn_locator_f[const] = np.array( [i,j] )
    #    i = j+1
    #for const in cbal: 
    #    j = i + eqn_locator_jac[const][1]-eqn_locator_jac[const][0]
    #    eqn_locator_f[const] = np.array( [i,j] ) 
    #    i = j+1
    #nf = i-1
    
    i = 1
    for const in f_list:
        j = i + eqn_locator_jac[const][1]-eqn_locator_jac[const][0]
        eqn_locator_f[const] = np.array( [i,j] )
        if i == j:
            idx2eqn[int(i)] = const
        else:
            for k in range(int(i),int(j)+1):
                idx2eqn[int(k)] = const + '[%i]'%(k+1-i)
        i = j+1
    nf = i-1
    
    print(idx2eqn)
    print(eqn_locator_f)
    
    
    #print(eqn_locator_f)
    #i = 1
    #for const in vle: 
    #    j = i + eqn_locator_jac[const][1]-eqn_locator_jac[const][0]
    #    eqn_locator_g[const] = np.array( [i,j] )
    #    i = j+1
    #for const in lflw:
    #    j = i + eqn_locator_jac[const][1]-eqn_locator_jac[const][0]
    #    eqn_locator_g[const] = np.array( [i,j] )
    #    i = j+1
    #for const in vflw:
    #    j = i + eqn_locator_jac[const][1]-eqn_locator_jac[const][0]
    #    eqn_locator_g[const] = np.array ( [i,j] ) 
    #    i = j+1
    #ng = i-1
    
    # Build index-finder for algebraic equations in g
    i = 1
    for const in g_list:
        j = i + eqn_locator_jac[const][1] - eqn_locator_jac[const][0]
        eqn_locator_g[const] = np.array( [i,j] )
        i = j+1 
    ng = i-1
    
    #i = 1
    #for var in Mvar:
    #    # i is incremented twice more to account for holdups at the 
    #    # reboiler and condenser, which are not included in the 
    #    # .col file
    #    i = i+1
    #    j = i + var_locator_jac[var][1]-var_locator_jac[var][0] 
    #    var_locator_z[var] = np.array( [i,j] )
    #    i = j+1
    #    i = i+1
    #for var in xvar: 
    #    j = i + var_locator_jac[var][1]-var_locator_jac[var][0]
    #    var_locator_z[var] = np.array( [i,j] )
    #    i = j+1
    #nz = i-1
    
    # Build index-finder for differential variables in z 
    i = 1
    for var in z_list:
        j = i + var_locator_jac[var][1] - var_locator_jac[var][0]
        var_locator_z[var] = np.array( [i,j] )
        if i == j:
            idx2state[int(i)] = var 
        else:
            for k in range(int(i),int(j)+1):
                idx2state[int(k)] = var + '[%i]'%(k+1-i)
        i = j+1
    nz = i-1

    print(var_locator_z)
    print(idx2state)
    
    #i = 1
    #for var in yvar:
    #    j = i + var_locator_jac[var][1]-var_locator_jac[var][0] 
    #    var_locator_y[var] = np.array( [i,j] )
    #    i = j+1
    #for var in Lvar: 
    #    j = i + var_locator_jac[var][1]-var_locator_jac[var][0]
    #    var_locator_y[var] = np.array( [i,j] )
    #    i = j+1
    #for var in Vvar: 
    #    j = i + var_locator_jac[var][1]-var_locator_jac[var][0]
    #    var_locator_y[var] = np.array( [i,j] )
    #    i = j+1
    #ny = i-1
    
    
    # Build index-finder for algebraic variables in y
    i = 1
    for var in y_list:
        j = i + var_locator_jac[var][1] - var_locator_jac[var][0]
        var_locator_y[var] = np.array( [i,j] )
        i = j+1
    ny = i-1
    
    # TODO: fix problem where stale variables are given a spot in 
    # y vector. This happens because I naively construct y_list
    # without checking for stale variables
    # ... could check if variable was assigned a non-zero value in the jacobian...
    # works because jacobian locator reads names directly from .col file
    # which does not include stale variables.
    
    #print('differential equations:')
    #for e in f_list: print(e,eqn_locator_f[e])
    #print('algebraic equations:')
    #for e in g_list: print(e,eqn_locator_g[e])
    #print('differential variables:')
    #for v in z_list: print(v,var_locator_z[v])
    #print('algebraic variables:')
    #for v in y_list: print(v,var_locator_y[v])
    
    
    nf = int(nf)
    nz = int(nz)
    ng = int(ng)
    ny = int(ny)
    print(nf,nz,ng,ny)
    
    
    dfdz = np.zeros((nf,nz))
    dfdy = np.zeros((nf,ny))
    dgdz = np.zeros((ng,nz))
    dgdy = np.zeros((ng,ny))
    
    # the following map equations to the variables they contain
    # 
    # ...
    # don't even end up using these, this method doesn't generalize well to different models
    # should be an easier way to do this with pyomo constraints/variables, if I think
    # its necessary
    # ...
    z_in_f = {}
    y_in_f = {}
    z_in_g = {}
    y_in_g = {}
    
    #for eqn in f:
    #    for # need list of variables contained in each equation
    #    dfdz[eqn_locator_f[eqn][0]:eqn_locator_f[eqn][1],] = \
    
    #for eqn in f:
    #    if eqn in mbal: 
    #        z_in_f[eqn] = [v for v in Mvar]
    #        y_in_f[eqn] = [v for v in Lvar + Vvar]
    #    if eqn in cbal:
    #        z_in_f[eqn] = [v for v in Mvar + xvar]
    #        y_in_f[eqn] = [v for v in yvar + Lvar + Vvar]
    #
    #for eqn in g: 
    #    if eqn in vle:
    #        z_in_g[eqn] = [v for v in xvar]
    #        y_in_g[eqn] = [v for v in yvar]
    #    if eqn in vflw:
    #        z_in_g[eqn] = []
    #        y_in_g[eqn] = [v for v in Vvar]
    #    if eqn in lflw:
    #        z_in_g[eqn] = [v for v in Mvar]
    #        y_in_g[eqn] = [v for v in Lvar]
    
    #print(z_in_f)
    #print(y_in_g) 
    
    '''
    #
    #for eqn in f:
    #    #print(eqn)
    #    for var in z_in_f[eqn]:
    #        #print(var) 
    #        dfdz[ int(eqn_locator_f[eqn][0])-1:int(eqn_locator_f[eqn][1]), \
    #                int(var_locator_z[var][0])-1:int(var_locator_z[var][1]) ] = \
    #        const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
    #                int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
    #    for var in y_in_f[eqn]:
    #        dfdy[ int(eqn_locator_f[eqn][0])-1:int(eqn_locator_f[eqn][1]), \
    #                int(var_locator_y[var][0])-1:int(var_locator_y[var][1]) ] = \
    #        const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
    #                int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
    #        
    #for eqn in g:
    #    #print(eqn)
    #    for var in z_in_g[eqn]:
    #        #print(var) 
    #        dgdz[ int(eqn_locator_g[eqn][0])-1:int(eqn_locator_g[eqn][1]), \
    #                int(var_locator_z[var][0])-1:int(var_locator_z[var][1]) ] = \
    #        const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
    #                int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
    #    for var in y_in_g[eqn]:
    #        dgdy[ int(eqn_locator_g[eqn][0])-1:int(eqn_locator_g[eqn][1]), \
    #                int(var_locator_y[var][0])-1:int(var_locator_y[var][1]) ] = \
    #        const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
    #                int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
    '''
    
    # constructing matrices
    # (probably) irrelevant but interesting question:
    # does the following code exploit locality well? 
    # intuition is no, because I grab data that can span
    # multiple matrix rows at once, and may revisit the same 
    # row many times. But is this worse 
    # than making more jumps between my source and 
    # target matrices? ... need to review ...
    for eqn in f_list:
        #print(eqn)
        for var in z_list:
            #print(var) 
            dfdz[ int(eqn_locator_f[eqn][0])-1:int(eqn_locator_f[eqn][1]), \
                    int(var_locator_z[var][0])-1:int(var_locator_z[var][1]) ] = \
            const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
                    int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
        for var in y_list:
            dfdy[ int(eqn_locator_f[eqn][0])-1:int(eqn_locator_f[eqn][1]), \
                    int(var_locator_y[var][0])-1:int(var_locator_y[var][1]) ] = \
            const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
                    int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
            
    for eqn in g_list:
        #print(eqn)
        for var in z_list:
            #print(var) 
            dgdz[ int(eqn_locator_g[eqn][0])-1:int(eqn_locator_g[eqn][1]), \
                    int(var_locator_z[var][0])-1:int(var_locator_z[var][1]) ] = \
            const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
                    int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
        for var in y_list:
            dgdy[ int(eqn_locator_g[eqn][0])-1:int(eqn_locator_g[eqn][1]), \
                    int(var_locator_y[var][0])-1:int(var_locator_y[var][1]) ] = \
            const_jac[ int(eqn_locator_jac[eqn][0])-1:int(eqn_locator_jac[eqn][1]), \
                    int(var_locator_jac[var][0])-1:int(var_locator_jac[var][1]) ]
    
    
    np.savetxt('dfdz.txt',dfdz,fmt = '%1.2f')
    np.savetxt('dfdy.txt',dfdy,fmt = '%1.2f')
    np.savetxt('dgdz.txt',dgdz,fmt = '%1.2f')
    np.savetxt('dgdy.txt',dgdy,fmt = '%1.2f')
    
    #print(dgdy)
    #np.savetxt('dgdy.txt',dgdy,fmt = '%1i')
    dgdy_inv = la.inv(dgdy)
    #dfdg = np.matmul(dfdy,dgdy_inv)
    
    dydz = np.matmul(dgdy_inv,dgdz)
    jacobian = dfdz - np.matmul(dfdy,dydz) 

    # compute idx2state dictionary to return as well
    # (index starting with 1 to match intuition about matrices)
    #for idx in range(1,nz+1):


    return jacobian, idx2state, idx2eqn
#print(np.matmul(dfdy,dydz))
#np.savetxt('djac_test',jacobian-dfdz,fmt = '%1.2f')


def main(): 
    jacobian, idx2state, idx2eqn = get_jacobian('jacobi_debug.in')
    nz = len(jacobian)
    lam, evec = la.eig(jacobian)
    N,s,LT = la.svd(jacobian)
    print(s)
    print(LT)
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
    
    
    print(x_1l)
    
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

if __name__ == '__main__':
    main()
