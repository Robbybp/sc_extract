from pyomo.environ import * 
from pyomo.dae import ContinuousSet, DerivativeVar
import numpy as np
from ss_mod import m_ss 
from ol_mod import m_ol

# initialize

#m.xE_IPA[1].set_value(0.0016)
#m.xE_IPA[2].set_value(0.0038)
#m.xE_IPA[3].set_value(0.0070)
#m.xE_IPA[4].set_value(0.0118)
#
#m.yE_IPA[1].set_value(0.000489)
#m.yE_IPA[2].set_value(0.001228)
#m.yE_IPA[3].set_value(0.002342)
#m.yE_IPA[4].set_value(0.004024)
#
#m.xE_F_IPA.set_value(0.019)
#m.yE_S_IPA.set_value(0.0)
#
#m.FFE.set_value(1.8)
#m.FSE.set_value(7.75)
#
#m.TE.set_value(313.0)
#
#m.KF_CO2.set_value(1.3)
#m.KF_IPA.set_value(0.0125)

ipopt = SolverFactory('ipopt')
ipopt.solve(m_ss,tee=True)

m.write('ss_mod.nl')#,io_options={symbolic_solver_labels:True})

with open('m_ss.txt', 'w') as f:
    m.display(ostream=f)
