from pyomo.environ import *
from pyomo.dae import ContinuousSet, DerivativeVar
import numpy as np
from ss_mod import m as m_ss

def solve_ss(m,u_vec,w_vec):
    # only supports V/L/Fcool control now, may be changed in the future
    m.U1_V.set_value(u_vec[1])
    m.U2_L.set_value(u_vec[2])
    m.U3_Fcool.set_value(u_vec[3])

    m.W1_MK.set_value(w_vec[1])
    m.W2_Tcool.set_value(w_vec[2])
    m.W3_F.set_value(w_vec[3])
    m.W4_xF.set_value(w_vec[4])

    ipopt = SolverFactory('ipopt')
    ipopt.solve(m,tee=True)

    return m

# temporary testing:
#
#u_vec = { 1 : 0.7185,
#          2 : 1.0,
#          3 : 22.661805 }
#
#w_vec = { 1 : 0.322074,
#          2 : 293,
#          3 : 1.8,
#          4 : 0.019 }
#
#m = solve_ss(m_ss,u_vec,w_vec)
#
#with open('m_ss.txt','w') as f:
#    m.display(ostream=f)
