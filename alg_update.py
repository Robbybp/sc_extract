from pyomo.environ import *
from pyomo.dae import ContinuousSet, DerivativeVar
import numpy as np
#from ol_mod import m as m_ol

# Function to update the algebraic states of model at time t
# when the differential states, inputs, and disturbances
# are already specified at this time.
def alg_update(m,t):

    # Extractor parameters 
    m.TC_sh_out[t].set_value(2*m.TC_sh_ave[t].value-m.TC_sh_in[t].value)
    m.TC_tb_out[t].set_value(2*m.TC_tb_ave[t].value-m.TC_tb_in[t].value)
    m.TE[t].set_value(m.TC_tb_out[t].value)

    m.KE_IPA[t].set_value(-0.002*m.TE[t].value+1.016)

    # flash inlet enthalpies
    H_CO2_l_298 = ((195.67-310*10.18)*298.15 + 10.18/2*298.15**2)/1000
    m.H_CO2_E[t].set_value(((195.67-310*10.18)*m.TE[t].value + 10.18/2*m.TE[t].value**2)/1000 - H_CO2_l_298)

    H_IPA_l_298 = (165.6-311.6*0.756)*298.15 + 0.756/2*298.15**2
    m.H_IPA_E[t].set_value(((165.6-311.6*0.756)*m.TE[t].value + 0.756/2*m.TE[t].value**2-H_IPA_l_298)/1000)

    # Extractor vapor liquid equilibrium
    m.yE_IPA[4,t].set_value((m.KE_IPA[t].value/m.EMVE.value)*(m.xE_IPA[4,t].value-m.xE_F_IPA[t].value*(1-m.EMVE.value)))
    for i in [3,2,1]:    
        m.yE_IPA[i,t].set_value((m.KE_IPA[t].value/m.EMVE.value)*(m.xE_IPA[i,t].value-m.xE_IPA[i+1,t].value*(1-m.EMVE.value)))

    # Flash VLE
    m.KF_CO2.set_value(1.3)
    m.KF_IPA.set_value(0.0125)
    m.xF_CO2[t].set_value( 1/(m.KF_CO2.value - m.KF_IPA.value) * (1 - m.KF_IPA.value))
    m.yF_CO2[t].set_value( m.KF_CO2.value * m.xF_CO2[t].value )
    # split fraction:
    m.qF[t].set_value( (1 - m.yE_IPA[m.NTE.value,t].value - m.xF_CO2[t].value)/(m.yF_CO2[t].value - m.xF_CO2[t].value) )
    
    # Flash temperature
    # might need to solve an implicit function (sub-problem) here
    m.TF[t].set_value(350)
    
    A = 25.0
    B = 55.19
    C = -33.69
    D = 7.948
    E = -0.1366
    F = -403.61
    H = -393.51
    m.H_CO2_v[t].set_value( A*(m.TF[t].value/1000) + B/2* (m.TF[t].value/1000)**2 + \
            C/3*(m.TF[t].value/1000)**3 + D/4*(m.TF[t].value)**4 - \
            E/(m.TF[t].value/1000) +F - H)
    m.H_CO2_l[t].set_value( ((195.69 - 310*10.18)*m.TF[t].value + \
            10.18/2*m.TF[t].value**2)/1000 - H_CO2_l_298 )
    H_IPA_v_298 = ((89.74 - 300*0.219)*298.15 + 0.219/2*298.15**2 )/1000
    m.H_IPA_v[t].set_value( ((89.74 - 300*0.219)*m.TF[t].value + 0.219/2*m.TF[t].value**2)/1000 \
            - H_IPA_v_298 )
    m.H_IPA_l[t].set_value( ((165.6-311.6*0.756)*m.TF[t].value + 0.756/2*m.TF[t].value**2 \
            + H_IPA_l_298)/1000 )






    return m

