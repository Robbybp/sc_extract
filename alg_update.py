from pyomo.environ import *
from pyomo.dae import ContinuousSet, DerivativeVar
import numpy as np
from math import log10
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



    # Flash VLE
    m.KF_CO2.set_value(m.KF_CO2_nom.value)
    m.KF_IPA.set_value(m.KF_IPA_nom.value)
    m.xF_CO2[t].set_value( 1/(m.KF_CO2.value - m.KF_IPA.value) * (1 - m.KF_IPA.value))
    m.yF_CO2[t].set_value( m.KF_CO2.value * m.xF_CO2[t].value )
    m.xF_IPA[t].set_value( 1-m.xF_CO2[t].value )
    m.yF_IPA[t].set_value( 1-m.yF_CO2[t].value )


    ###########################################
    # Block of code to calculate qF, yE_S_IPA # 
    ###########################################

# This code was intended to solve for qF,yE_S_IPA in the original model m
# It is unused due to an error adding constraints (probably because model is already discretized)

#    for const in m.component_objects(Constraint):
#        const.deactivate()
#    for var in m.component_objects(Var):
#        print(var)
#        #if 'fix' in dir(var):
#        if not isinstance(var,DerivativeVar):
#            var.fix()
#
#    m.qF.unfix()
#    m.yE_S_IPA.unfix()


# sub-model 1
    submodel_1 = ConcreteModel()
# variables, solve for:
    # vapor fraction after flash:
    submodel_1.qF = Var(initialize=1)
    # Solvent IPA fraction, fed to extractor:
    submodel_1.yS = Var(initialize=0)
# parameters (some of which were variables in the full open-loop model)
    # vapor and liquid CO2 fractions after the flash
    submodel_1.yF  = Param(initialize=m.yF_CO2[t].value)
    submodel_1.xF  = Param(initialize=m.xF_CO2[t].value)
    # extractor IPA VLE coefficient
    submodel_1.K   = Param(initialize=m.KE_IPA[t].value)
    # extractor tray efficiency 
    submodel_1.E   = Param(initialize=m.EMVE)
    # liquid IPA composition in the extractor
    submodel_1.x1  = Param(initialize=m.xE_IPA[1,t].value)
    submodel_1.x2  = Param(initialize=m.xE_IPA[2,t].value)
    submodel_1.x3  = Param(initialize=m.xE_IPA[3,t].value)
    submodel_1.x4  = Param(initialize=m.xE_IPA[4,t].value)
    # boil-up flow rate in stripper
    submodel_1.V   = Param(initialize=m.FV_S[t].value)
    # make up flow rate in condenser
    submodel_1.FMK = Param(initialize=m.FMK[t].value)
    # CO2 composition in the top (4th) tray of the stripper
    submodel_1.yS4 = Param(initialize=m.yS_CO2[4,t].value)
    # Total flow rate of the feed to the stripper
    submodel_1.FS  = Param(initialize=m.FF_S[t].value)

# Equations to solve
    def sm1_const_1_rule(m):
        return m.qF == 1/(m.yF-m.xF)*(1 - m.xF -\
                (m.K*m.E*m.x4 + m.K*m.E*m.x3*(1-m.E) + m.K*m.E*m.x2*(1-m.E)**2 + \
                m.K*m.E*m.x1*(1-m.E)**3 + m.yS*(1-m.E)**4 ) ) 
    submodel_1.const_1 = Constraint(rule=sm1_const_1_rule)

    def sm1_const_2_rule(m):
        return 1-m.yS == ( (m.V + m.qF*m.FS)*m.yS4 + m.FMK )/ \
                         ( m.FMK + m.V + m.qF*m.FS )
    submodel_1.const_2 = Constraint(rule=sm1_const_2_rule)

    with open('submodel_1_init.txt','w') as f:
        submodel_1.display(ostream=f)

    submodel_1.write('sub1.nl')

    ipopt = SolverFactory('ipopt') 
    ipopt.solve(submodel_1,tee=True)

    with open('submodel_1.txt','w') as f:
        submodel_1.display(ostream=f)


# Adding this constraint to the original model caused an error
    #def yS_const_1_rule(m):
    #    return m.qF[0] == (1 - m.xF_CO2[0])/(m.yF_CO2[0] - m.xF_CO2[0]) - 1/(m.xF_CO2[0]-m.yF_CO2[0])*\
    #            m.EMVE*m.KE_IPA[0]*m.xE_IPA[4,0] + \
    #            m.EMVE*m.KE_IPA[0]*m.xE_IPA[3,0]*(1-m.EMVE) + \
    #            m.EMVE*m.KE_IPA[0]*m.xE_IPA[2,0]*(1-m.EMVE)**2 + \
    #            m.EMVE*m.KE_IPA[0]*m.xE_IPA[1,0]*(1-m.EMVE)**3 + \
    #            m.EMVE*m.KE_IPA[0]*m.yE_S_IPA[0]*(1-m.EMVE)**4
    #m.yS_const_1 = Constraint(rule=yS_const_1_rule)


# Code to unfix variables and re-activate constraints, again, unused
#    for const in m.component_objects(Constraint):
#        const.activate()
#    #m.yS_const_1.deactivate()
#    for var in m.component_objects(Var):
#        if not isinstance(var,DerivativeVar):
#            var.unfix()

    m.yE_S_IPA[t].set_value(submodel_1.yS.value)
    m.qF[t].set_value(submodel_1.qF.value)

#############################################################################
    
    # Calculate y_stripper from x_stripper
    m.KS_IPA.set_value(m.KS_IPA_nom.value)
    m.yS_CO2[1,t].set_value( 1- m.KS_IPA.value*(1 - m.xS_CO2[1,t].value) )
    m.yS_CO2[2,t].set_value( 1- m.KS_IPA.value*m.EMVS.value * (1 - m.xS_CO2[2,t].value) - \
            (1 - m.yS_CO2[1,t].value) * (1 - m.EMVS.value) )
    m.yS_CO2[3,t].set_value( 1- m.KS_IPA.value*m.EMVS.value * (1 - m.xS_CO2[3,t].value) - \
            (1 - m.yS_CO2[2,t].value) * (1 - m.EMVS.value) )
    m.yS_CO2[4,t].set_value( 1- m.KS_IPA.value*m.EMVS.value * (1 - m.xS_CO2[4,t].value) - \
            (1 - m.yS_CO2[3,t].value) * (1 - m.EMVS.value) )

    m.xS_CO2[5,t].set_value( 1 - m.yE_S_IPA[t].value )
    m.yS_CO2[5,t].set_value( m.xS_CO2[5,t].value )


    # Extractor VLE
    m.yE_IPA[1,t].set_value( m.KE_IPA[t].value*m.EMVE*m.xE_IPA[1,t].value + \
            m.yE_S_IPA[t].value*(1-m.EMVE) )
    m.yE_IPA[2,t].set_value( m.KE_IPA[t].value*m.EMVE*m.xE_IPA[2,t].value + \
            m.yE_IPA[1,t].value*(1-m.EMVE) )
    m.yE_IPA[3,t].set_value( m.KE_IPA[t].value*m.EMVE*m.xE_IPA[3,t].value + \
            m.yE_IPA[2,t].value*(1-m.EMVE) )
    m.yE_IPA[4,t].set_value( m.KE_IPA[t].value*m.EMVE*m.xE_IPA[4,t].value + \
            m.yE_IPA[3,t].value*(1-m.EMVE) )

    # flash inlet enthalpies
    H_CO2_l_298 = ((195.67-310*10.18)*298.15 + 10.18/2*298.15**2)/1000
    m.H_CO2_E[t].set_value(((195.67-310*10.18)*m.TE[t].value + 10.18/2*m.TE[t].value**2)/1000 - H_CO2_l_298)

    H_IPA_l_298 = (165.6-311.6*0.756)*298.15 + 0.756/2*298.15**2
    m.H_IPA_E[t].set_value(((165.6-311.6*0.756)*m.TE[t].value + 0.756/2*m.TE[t].value**2-H_IPA_l_298)/1000)

# Initialize flash temperature to temperature fed to the flash
    m.TF[t].set_value(m.TE[t].value)

    A = 25.0
    B = 55.19
    C = -33.69
    D = 7.948
    E = -0.1366
    F = -403.61
    H = -393.51
    m.H_CO2_v[t].set_value( A*(m.TF[t].value/1000) + B/2* (m.TF[t].value/1000)**2 + \
            C/3*(m.TF[t].value/1000)**3 + D/4*(m.TF[t].value/1000)**4 - \
            E/(m.TF[t].value/1000) +F - H)
    m.H_CO2_l[t].set_value( ((195.69 - 310*10.18)*m.TF[t].value + \
            10.18/2*m.TF[t].value**2)/1000 - H_CO2_l_298 )
    H_IPA_v_298 = ((89.74 - 300*0.219)*298.15 + 0.219/2*298.15**2 )/1000
    m.H_IPA_v[t].set_value( ((89.74 - 300*0.219)*m.TF[t].value + 0.219/2*m.TF[t].value**2)/1000 \
            - H_IPA_v_298 )
    m.H_IPA_l[t].set_value( ((165.6-311.6*0.756)*m.TF[t].value + 0.756/2*m.TF[t].value**2 \
            - H_IPA_l_298)/1000 )


    # Flash temperature
    # might need to solve an implicit function (sub-problem) here
########################################################################
    
# Submodel for solving for flash temperature
    submodel_2 = ConcreteModel()
# In this problem, variables are flash temperature and enthalpies 
    submodel_2.TF = Var(initialize=m.TE[t].value)
    submodel_2.H_CO2_l = Var(initialize=m.H_CO2_l[t].value)
    submodel_2.H_CO2_v = Var(initialize=m.H_CO2_v[t].value)
    submodel_2.H_IPA_l = Var(initialize=m.H_IPA_l[t].value)
    submodel_2.H_IPA_v = Var(initialize=m.H_IPA_v[t].value)
# Parameters ported from the main model:
    submodel_2.qF = Param(initialize=m.qF[t].value)
    # CO2 compositions in the liquid and vapor after the flash
    submodel_2.xF = Param(initialize=m.xF_CO2[t].value)
    submodel_2.yF = Param(initialize=m.yF_CO2[t].value)
    # IPA composition of supercritical stream entering flash from extractor
    submodel_2.yE = Param(initialize=m.yE_IPA[4,t].value)
    # Enthalpies of each component entering flash
    submodel_2.H_CO2_E = Param(initialize=m.H_CO2_E[t].value)
    submodel_2.H_IPA_E = Param(initialize=m.H_IPA_E[t].value)
# Equations to satisfy in this sub-model
    # four enthalpy definition equations:
    def sm2_const_1_rule(m):
        return m.H_CO2_v == \
            A*(m.TF/1000) + B/2* (m.TF/1000)**2 + \
            C/3*(m.TF/1000)**3 + D/4*(m.TF/1000)**4 - \
            E/(m.TF/1000) +F - H
    submodel_2.const_1 = Constraint(rule=sm2_const_1_rule)

    def sm2_const_2_rule(m):
        return m.H_CO2_l == \
            ((195.69 - 310*10.18)*m.TF + \
            10.18/2*m.TF**2)/1000 - H_CO2_l_298 
    submodel_2.const_2 = Constraint(rule=sm2_const_2_rule)

    H_IPA_v_298 = ((89.74 - 300*0.219)*298.15 + 0.219/2*298.15**2 )/1000
    def sm2_const_3_rule(m):
        return m.H_IPA_l == \
            ((165.6-311.6*0.756)*m.TF + 0.756/2*m.TF**2 \
            - H_IPA_l_298)/1000 
    submodel_2.const_3 = Constraint(rule=sm2_const_3_rule)

    def sm2_const_4_rule(m):
        return m.H_IPA_v == \
            ((89.74 - 300*0.219)*m.TF + 0.219/2*m.TF**2)/1000 \
            - H_IPA_v_298 
    submodel_2.const_4 = Constraint(rule=sm2_const_4_rule)

    # one enthalpy balance equation:
    def sm2_const_5_rule(m):
        return m.yE*m.H_IPA_E + (1-m.yE)*m.H_CO2_E == \
                (1-m.qF)*(m.xF*m.H_CO2_l + (1-m.xF)*m.H_IPA_l) +\
                m.qF*(m.yF*m.H_CO2_v + (1-m.yF)*m.H_CO2_v)
    submodel_2.const_5 = Constraint(rule=sm2_const_5_rule)

    with open('submodel_2_init.txt','w') as f:
        submodel_2.display(ostream=f)

    ipopt.solve(submodel_2,tee=True)

    with open('submodel_2.txt','w') as f:
        submodel_2.display(ostream=f)

# Port results back into base model:
    
    m.H_CO2_l[t].set_value( submodel_2.H_CO2_l.value )
    m.H_CO2_v[t].set_value( submodel_2.H_CO2_v.value )
    m.H_IPA_l[t].set_value( submodel_2.H_IPA_l.value )
    m.H_IPA_v[t].set_value( submodel_2.H_IPA_v.value )

    m.TF[t].set_value( submodel_2.TF.value )

##################################################################    

    

    m.P1_F.set_value( m.P1_F_nom.value )
    m.P2_F.set_value( m.P2_F_nom.value )
    m.PS.set_value( m.P2_F.value )
    # nominal value of disturbance used because optimization problem will be solved with nominal
    # values in mind
    # ... but this is an open loop problem... should real values of disturbance be used...
    # I think so... and then they should be changd back to nominal values for MPC 
    m.FSE[t].set_value( 1/(1-m.qF[t].value)*(m.W1_MK[t].value + m.U1_V[t].value - m.U2_L[t].value))
    m.FFF[t].set_value( m.FSE[t].value )

    # Stripping Column flow rates:
    m.FF_S[t].set_value( m.FFF[t].value )
    m.FT_S[t].set_value( m.FV_S[t].value + m.qF[t].value*m.FF_S[t].value )
    m.FB_S[t].set_value( m.FL_S[t].value + (1-m.qF[t].value)*m.FF_S[t].value - m.FV_S[t].value )
    m.F_comp[t].set_value( m.FT_S[t].value + m.FMK[t].value )
    m.FD_S[t].set_value( m.F_comp[t].value - m.FL_S[t].value )

    m.xF_S[t].set_value( m.qF[t].value*m.yF_CO2[t].value + (1-m.qF[t].value)*m.xF_CO2[t].value )
    
    # Condenser variables
    m.P1_cond.set_value( m.P2_F.value )
    m.P2_cond.set_value( m.P2_cond_nom.value )
    m.T1_cond[t].set_value( m.TF[t].value )
    m.V1_cond[t].set_value( (1/m.P1_cond.value)*m.F_comp[t].value*0.08206*m.T1_cond[t].value )
    m.V2_cond[t].set_value( (m.P1_cond.value*m.V1_cond[t].value**1.4/m.P2_cond.value)**(1/1.4) )
    m.T2_cond[t].set_value( m.P2_cond.value*m.V2_cond[t].value/(m.F_comp[t].value*0.08206) )

    # Reboiler variables
    m.P_B.set_value( m.PS.value )
    m.Pvap_B_IPA[t].set_value( m.P_B.value*(1-m.yS_CO2[1,t].value)/(1-m.xS_CO2[1,t].value) )
    m.TB_sh[t].set_value( 273.15 + (1360 + 197.6*log10(m.Pvap_B_IPA[t].value*101.325) - \
            6.866*197.6)/(6.866 - log10(m.Pvap_B_IPA[t].value*101.325) ) )
    m.dH_IPA.set_value( m.dH_IPA_nom.value )
    m.dH_CO2.set_value( m.dH_CO2_nom.value )
    m.dH_B[t].set_value(m.yS_CO2[1,t].value*m.dH_CO2.value + (1-m.yS_CO2[1,t].value)*m.dH_IPA.value)
    m.CpB_CO2[t].set_value( m.CpB_CO2_nom.value )
    m.TB_tb_in[t].set_value( m.T2_cond[t].value )
    # TB_tb_ave known, since it's a differential variable
    m.TB_tb_out[t].set_value( 2*m.TB_tb_ave[t].value - m.TB_tb_in[t].value )
    m.F_B[t].set_value( m.FV_S[t].value*m.dH_B[t].value/\
            ((m.TB_tb_in[t].value - m.TB_tb_out[t].value)*m.CpB_CO2[t].value) )
    m.F_B_bypass[t].set_value( m.FD_S[t].value - m.F_B[t].value )

    # Cooler variables
    # TC_sh_ave and TC_tb_ave differential, so we may assume access to them in this function

    # should be the same as FSE calculated above:
    m.F_C[t].set_value( m.F_B[t].value + m.F_B_bypass[t].value )
    m.TC_tb_in[t].set_value( 1/m.F_C[t].value*(\
            m.TB_tb_out[t].value*m.F_B[t].value + m.T2_cond[t].value*m.F_B_bypass[t].value ) )
    m.TC_tb_out[t].set_value( 2*m.TC_tb_ave[t].value - m.TC_tb_in[t].value )
    m.TC_sh_out[t].set_value( 2*m.TC_sh_ave[t].value - m.TC_sh_in[t].value )
    m.CpC_CO2[t].set_value( m.CpB_CO2[t].value )
    m.CpC_H2O[t].set_value( 36.9 )
    m.yC_CO2[t].set_value( m.xS_CO2[m.NTS,t].value )
    m.yC_IPA[t].set_value( m.yE_S_IPA[t].value )


     
    





    return m
