from pyomo.environ import * 
from pyomo.dae import ContinuousSet, DerivativeVar
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
    m.xE_IPA_0[i]=(m_ss.xE_IPA[i].value)
for i in m.set_S_dyn:
    m.xS_CO2_0[i]=(m_ss.xS_CO2[i].value)
m.TB_tb_0=(m_ss.TB_tb_ave.value)
m.TC_tb_0=(m_ss.TC_tb_ave.value)
m.TC_sh_0=(m_ss.TC_sh_ave.value)

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
for t in m.t:
    m.U1_V[t] = m_ss.U1_V.value
    m.U2_L[t] = m_ss.U2_L.value
    m.U3_Fcool[t] = m_ss.U3_Fcool.value

    m.W1_MK[t] = m_ss.W1_MK.value
    m.W2_Tcool[t] = m_ss.W2_Tcool.value
    m.W3_F[t] = m_ss.W3_F.value
    m.W4_xF[t] = m_ss.W4_xF.value
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

m = alg_update(m,0)

with open('m_ol_init.txt','w') as f:
    m.display(ostream=f)

ipopt = SolverFactory('ipopt')
ipopt.solve(m,tee=True)

#m_ol.write('ol_mod.nl')#,io_options={symbolic_solver_labels:True})


with open('m_ol.txt', 'w') as f:
    m.display(ostream=f)
