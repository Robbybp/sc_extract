from pyomo.environ import * 
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.core.base.sets import _SetProduct
import numpy as np
from ss_mod import m as m_ss 
from ol_mod import m
from ss_fcn import solve_ss
from alg_update import alg_update

# Solve steady state model

u_vec = { 1 : 0.7185,
          2 : 1.0,
          3 : 22.661805 }

w_vec = { 1 : 0.322074,
          2 : 293,
          3 : 1.8,
          4 : 0.019 }

m_ss = solve_ss(m_ss,u_vec,w_vec)


with open('m_ss.txt','w') as f:
    m_ss.display(ostream=f)

m_ss.write('ss_mod.nl')

# Set differential state IC parameters

for i in m.set_E_Tray:
    m.xE_IPA_0[i].set_value(m_ss.xE_IPA[i].value)
for i in m.set_S_dyn:
    m.xS_CO2_0[i].set_value(m_ss.xS_CO2[i].value)
m.TB_tb_0.set_value(m_ss.TB_tb_ave.value)
m.TC_tb_0.set_value(m_ss.TC_tb_ave.value)
m.TC_sh_0.set_value(m_ss.TC_sh_ave.value)

# Initialize differential states 

for i in m.set_E_Tray:
    m.xE_IPA[i,0].set_value(m.xE_IPA_0[i].value)
for i in m.set_S_dyn:
    m.xS_CO2[i,0].set_value(m.xS_CO2_0[i].value)
m.TB_tb_ave[0].set_value(m.TB_tb_0.value)
m.TC_tb_ave[0].set_value(m.TC_tb_0.value)
m.TC_sh_ave[0].set_value(m.TC_sh_0.value)

# Initialize input and disturbance values in dynamic model

# this is is where the values used in the open loop simulation are set
# To run with time-varying inputs or disturbances,
# change these parameters to be time-varying
#print(m.t.get_finite_elements())
for t in m.t:
    m.U1_V[t].set_value(m_ss.U1_V.value)
    m.U2_L[t].set_value(m_ss.U2_L.value)
    m.U3_Fcool[t].set_value(m_ss.U3_Fcool.value)

    m.W1_MK[t].set_value(m_ss.W1_MK.value)
    m.W2_Tcool[t].set_value(m_ss.W2_Tcool.value)
    m.W3_F[t].set_value(m_ss.W3_F.value)
    m.W4_xF[t].set_value(m_ss.W4_xF.value)
# now initialize input/disturbance variable values    
    m.FV_S[t].set_value(m.U1_V[t].value)
    m.FL_S[t].set_value(m.U2_L[t].value)
    m.F_C_H2O[t].set_value(m.U3_Fcool[t].value)

    m.FMK[t].set_value(m.W1_MK[t].value)
    m.TC_sh_in[t].set_value(m.W2_Tcool[t].value)
    m.FFE[t].set_value(m.W3_F[t].value)
    m.xE_F_IPA[t].set_value(m.W4_xF[t].value)

# Initialize algebraic states
# write function for this...

#m = alg_update(m,0)

#ol_xE = getattr(m,'xE_IPA')
#print(dir(ol_xE))
#print(dir(ol_xE.index_set()))
#print([v for v in ol_xE.index_set()])
#print(ol_xE.index_set().set_tuple)
#print(m.t in ol_xE.index_set().set_tuple)

#print([v for _, v in ol_xE.index_set()])
#print(m.t in ol_xE.index_set())
#ol_K = getattr(m,'KF_CO2')
#ol_F = getattr(m,'FFE')
#print(ol_F)
#print(ol_F.index_set())
#print(type(ol_F.index_set()))
#print(ol_K)
#print(ol_K.index_set())
#print([v for v in ol_K.index_set()])
#print([t for t in m.t])
# Copy steady state initialization for all time

#'''
for var_ss in m_ss.component_objects(Var, active=True):
#    print(var_ss)
#    print(type(var_ss))
#    print(str(var_ss))
#    print(type(str(var_ss)))
    #print(var_ss.value)
    var_name=str(var_ss)
    var_ol = getattr(m,var_name)
    if type(var_ss.index_set()) is _SetProduct:
        ss_index_sets = var_ss.index_set().set_tuple
    else:
        ss_index_sets = var_ss.index_set()
#    print(ss_index_sets)
    if type(var_ol.index_set()) is _SetProduct:
        ol_index_sets = var_ol.index_set().set_tuple
    else: 
        ol_index_sets = var_ol.index_set() 

#    print(ol_index_sets)
#    print('t in set? ',m.t in ol_index_sets)
#    print('t is idx? ',m.t == ol_index_sets)

#    print('var_ol: \t',var_ol)
    for index in var_ss:
#        print(index)
        var_ss_value = var_ss[index].value
        index_type = type(index)
        if index is None:
            if m.t in ol_index_sets:

                for t in m.t:
                    ol_index = t
                    var_ol[ol_index].set_value(var_ss_value)
#                print(ol_index)
            else: 
                var_ol.set_value(var_ss_value)
            continue 
        if index_type is int:
            for t in m.t:
                ol_index = (index,t)
                var_ol[ol_index].set_value(var_ss_value)
#            print(ol_index)
            continue
        if index_type is tuple:
            for t in m.t:
                ol_index = index + (t) 
                var_ol[ol_index].set_value(var_ss_value)
#            print(ol_index)
            continue

#        print(type(index))
#        print(var_ss[index])

#'''
#print(getattr(m_ss,'xE_IPA'))

#for index in var_ol:
#    print(index)
#    print(type(index))
#
#tuple1 = (1,2)
#print(tuple(tuple1))
#three = tuple(3)
#tuple2 = (tuple1,3)
#print(tuple2)

m = alg_update(m,0)

with open('m_ol_init.txt','w') as f:
    m.display(ostream=f)

ipopt = SolverFactory('ipopt')
ipopt.solve(m,tee=True)

m.write('ol_mod.nl')#,io_options={symbolic_solver_labels:True})


with open('m_ol.txt', 'w') as f:
    m.display(ostream=f)
