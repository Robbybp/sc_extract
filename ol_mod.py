from pyomo.environ import * 
from pyomo.dae import ContinuousSet, DerivativeVar
#from pyomo.dae import *
import numpy as np

m = ConcreteModel()

# Control options

# Model options
const_Cp_switch = 0

# Discretization options

m.t = ContinuousSet(bounds=(0,5)) # solve nmpc opt problem over 5 hours

# Parameters

# heat capacities at 40 C, 102 atm
m.Cp_H20_L = 4.186

# Extractor Parameters:
m.NTE = Param(initialize=4) # (integer)
m.ME = Param(initialize=0.25) # kmol, steady state holdup
# ^will probably need to rename these holdup parameters
m.EMVE = Param(initialize=0.75) # Murphy stage efficiency


# Stripper Parameters:
m.NTS = Param(initialize=5) # Number of trays, including reboiler and condenser  
m.NTFS = Param(initialize=3) # Feed tray; trays count from bottom up
m.MS = Param(initialize=0.05) # kmol, tray holdup
m.MBS = Param(initialize=0.25) # reboiler holdup
m.EMVS = Param(initialize=0.80) # 

# Reboiler Parameters:
m.LTB = Param(initialize=0.365) # tube length, meters
m.NTB = Param(initialize=34) # number of tubes
DTB = 0.5*0.0254
m.DTB = Param(initialize=DTB) # tube diameter, (inches->) meters
A_cr_B = m.NTB*0.25*3.1416*m.DTB**2
m.A_cr_B = Param(initialize=A_cr_B)
A_B = m.LTB*m.NTB*3.1416*m.DTB
m.A_B = Param(initialize=A_B)
Vol_B = m.LTB*m.A_cr_B
m.Vol_B = Param(initialize=Vol_B)
m.U_ov_B = Param(initialize=360) # overall heat transfer coefficient in reboiler, kJ/h/m^2/K

# Trim-cooler Parameters:
m.LTC = Param(initialize=1.15) # tube length, meters
m.NTC = Param(initialize=38) # number of tubes
DTC = 0.5*0.0254
m.DTC = Param(initialize=DTC) # tube diameter, (inches->) meters
ACrC = m.NTC*0.25*3.1416*m.DTC**2
m.ACrC = Param(initialize=ACrC)
A_C = m.NTC*3.1416*m.LTC*m.DTC # total area of heat transfer in cooler
m.A_C = Param(initialize=A_C)
m.A_C.pprint()
VolC = m.ACrC*m.LTC
m.VolC = Param(initialize=VolC) 
m.U_ov_C = Param(initialize=360) 

############
### Sets ###
############

m.set_E_Tray = Set(initialize=range(1,m.NTE.value+1))
m.set_E_int = Set(initialize=range(2,m.NTE.value))
m.set_E_2_NT = Set(initialize=range(2,m.NTE.value+1))
m.set_S_Tray = Set(initialize=range(1,m.NTS.value+1))
m.set_S_eq = Set(initialize=range(2,m.NTS.value))
m.set_S_dyn = Set(initialize=range(1,m.NTS.value))

####################################
### Initial/operating conditions ###
####################################

# Disturbances/parameters
m.FFE_nom = Param(initialize=1.8)
m.xE_F_nom = Param(initialize=0.019,mutable=True)
# Process variables
m.FSE_init = Param(initialize=7.75)
# ^won't be necessary when full model is connected
m.yE_S_sp = Param(initialize=0.0,mutable=True)
m.TE_sp = Param(initialize=313.0,mutable=True)

# Stripping Column
m.FMK_nom = Param(initialize=0.31875+0.00332521)
m.FVS_nom = Param(initialize=0.7185)
m.FLS_nom = Param(initialize=1.0)
m.FSS_nom = Param(initialize=7.75)
m.FTS_nom = Param(initialize=8.43)
m.F_comp_nom = Param(initialize=8.70975)
m.FDS_nom = Param(initialize=7.70975)
m.KS_IPA_nom = Param(initialize=0.0213)

### State initial conditions ###

m.xE_IPA_0 = Param(m.set_E_Tray,mutable=True,initialize={
    1 : 0.00214759,
    2 : 0.00419620,
    3 : 0.00728832,
    4 : 0.01195549})
m.xS_CO2_0 = Param(m.set_S_dyn,mutable=True,initialize={
    1 : 0.90508156,
    2 : 0.96946326,
    3 : 0.97021213,
    4 : 0.99574316})
m.TB_tb_0 = Param(mutable=True,initialize=413.4)
m.TC_tb_0 = Param(mutable=True,initialize=368.8)
m.TC_sh_0 = Param(mutable=True,initialize=313.7)

########################################
### Input and disturbance parameters ###
########################################

m.U1_V = Param(m.t,initialize=0.7185,default=0.7185,mutable=True)
m.U2_L = Param(m.t,initialize=1.0,mutable=True,default=1.0)
m.U3_Fcool = Param(m.t,initialize=22.661805,mutable=True,default=22.661805)
# unused for now, but may be used in the future:
m.U4_B = Param(m.t,initialize=0.322074,mutable=True,default=0.322074)
m.U5_D = Param(m.t,initialize=7.75,mutable=True,default=7.75)

#m.W1_S = Param(initialize=7.75,mutable=True)
#^Solvent flow rate not a disturbance, more of an input (as distillate from stripper)
m.W1_MK = Param(m.t,initialize=0.322074,mutable=True,default=0.322074)
m.W2_Tcool = Param(m.t,initialize=293,mutable=True,default=293)
m.W3_F = Param(m.t,initialize=1.8,mutable=True,default=1.8)
m.W4_xF = Param(m.t,initialize=0.019,mutable=True,default=0.019)

# nominal values
m.U1_V_nom = Param(m.t,initialize=0.7185,mutable=True,default=0.7185)
m.U2_L_nom = Param(m.t,initialize=1.0,mutable=True,default=1.0)
m.U3_Fcool_nom = Param(m.t,initialize=22.661805,mutable=True,default=22.661805)

#m.W1_S_nom = Param(m.t,initialize=7.75,mutable=True)
m.W1_MK_nom = Param(m.t,initialize=0.322074,mutable=True,default=0.322074)
m.W2_Tcool_nom = Param(m.t,initialize=293,mutable=True,default=293)
m.W3_F_nom = Param(m.t,initialize=1.8,mutable=True,default=1.8)
m.W4_xF_nom = Param(m.t,initialize=0.019,mutable=True,default=0.019)


#######################
### State variables ###
#######################

# Extractor:
# compositions below are that of IPA 
m.xE_IPA = Var(m.set_E_Tray,m.t,within=NonNegativeReals)
m.yE_IPA = Var(m.set_E_Tray,m.t,within=NonNegativeReals)
m.xE_F_IPA = Var(m.t,within=NonNegativeReals,initialize=0.019)
m.yE_S_IPA = Var(m.t,within=NonNegativeReals,initialize=0)
m.FFE = Var(m.t,within=NonNegativeReals,initialize=1.8) # Feed Flow into Extractor
# ^is this a parameter, state variable, or control variable? 
# A: uncertain parameter/disturbance
m.FSE = Var(m.t,within=NonNegativeReals,initialize=m.FSE_init) # Solvent Flow into Extractor
m.KE_IPA = Var(m.t,within=NonNegativeReals)
m.TE = Var(m.t,within=NonNegativeReals,initialize=313)

m.xE_IPAdot = DerivativeVar(m.xE_IPA,wrt=m.t,initialize=0)


# Flash:
m.FFF = Var(m.t,within=NonNegativeReals,initialize=m.FSE_init) # Flow of Feed to Flash
m.qF = Var(m.t,within=NonNegativeReals,initialize=0.999) # Vapor split fraction after flash
m.xF_IPA = Var(m.t,within=NonNegativeReals,initialize=0.2)
m.xF_CO2 = Var(m.t,within=NonNegativeReals)
m.yF_IPA = Var(m.t,within=NonNegativeReals,initialize=0.004)
m.yF_CO2 = Var(m.t,within=NonNegativeReals)
m.KF_CO2 = Var(within=NonNegativeReals)
m.KF_IPA = Var(within=NonNegativeReals)
m.HFF = Var(m.t,within=NonNegativeReals) # Enthalpy in the feed of the flash
m.HVF = Var(m.t,within=NonNegativeReals) # Enthalpy of the vapor after the flash
m.HLF = Var(m.t,within=NonNegativeReals) # Enthalpy of the liquid after the flash
m.Cp_CO2_l = Var(m.t,within=NonNegativeReals,initialize=226.21)
m.Cp_CO2_v = Var(m.t,within=NonNegativeReals,initialize=37.82)
m.Cp_IPA_l = Var(m.t,within=NonNegativeReals,initialize=166.66)
m.Cp_IPA_v = Var(m.t,within=NonNegativeReals,initialize=92.59)
m.H_CO2_l = Var(m.t,within=Reals,initialize=3.22) # why are the liquid streams at higher enthalpy?
m.H_CO2_v = Var(m.t,within=Reals,initialize=0.55)
m.H_CO2_E = Var(m.t,within=Reals,initialize=3.0)
m.H_IPA_l = Var(m.t,within=Reals,initialize=2.383)
m.H_IPA_v = Var(m.t,within=Reals,initialize=1.34)
m.H_IPA_E = Var(m.t,within=Reals,initialize=2.0)
m.TF = Var(m.t,within=NonNegativeReals,initialize=310)
m.P1_F = Var(within=NonNegativeReals,initialize=100)
m.P2_F = Var(within=NonNegativeReals,initialize=40)

# Stripping Column:
m.xS_CO2 = Var(m.set_S_Tray,m.t,within=NonNegativeReals,initialize=0.9)
m.yS_CO2 = Var(m.set_S_Tray,m.t,within=NonNegativeReals,initialize=0.9)
m.FF_S = Var(m.t,within=NonNegativeReals,initialize=m.FSE_init) # feed flow rate
m.FL_S = Var(m.t,within=NonNegativeReals,initialize=m.FLS_nom.value) # liquid reflux flow rate
m.FV_S = Var(m.t,within=NonNegativeReals,initialize=m.FVS_nom.value) # vapor boilup flow rate
m.FB_S = Var(m.t,within=NonNegativeReals) # bottoms product flow rate 
m.FMK = Var(m.t,within=NonNegativeReals,initialize=m.FMK_nom.value) # make up CO2 flow rate
#m.FS_S = Var(within=NonNegativeReals,initialize=m.FSS_nom.value) 
m.FT_S = Var(m.t,within=NonNegativeReals,initialize=m.FTS_nom.value) # flow off top of column
m.FD_S = Var(m.t,within=NonNegativeReals,initialize=m.FDS_nom.value) # distillate (solvent) flow rate
m.F_comp = Var(m.t,within=NonNegativeReals,initialize=m.F_comp_nom.value) # flow rate through compresser 
m.xF_S = Var(m.t,within=NonNegativeReals) # CO2 fraction of combined (liq & vap) feed to stripper
m.KS_CO2 = Var(within=NonNegativeReals)
m.KS_IPA = Var(within=NonNegativeReals,initialize=m.KS_IPA_nom.value)
m.TS = Var(m.t,within=NonNegativeReals)
m.PS = Var(within=NonNegativeReals)

m.xS_CO2dot = DerivativeVar(m.xS_CO2,wrt=m.t,initialize=0)

# Condenser:
m.P1_cond = Var(within=NonNegativeReals)
m.P2_cond = Var(within=NonNegativeReals)
m.T1_cond = Var(m.t,within=NonNegativeReals)
m.T2_cond = Var(m.t,within=NonNegativeReals)
m.V1_cond = Var(m.t,within=NonNegativeReals)
m.V2_cond = Var(m.t,within=NonNegativeReals)

# Reboiler:
m.rhoB_CO2 = Var(m.t,within=NonNegativeReals) # Density in kg/m^3
m.CpB_CO2 = Var(m.t,within=NonNegativeReals)
m.TB_sh = Var(m.t,within=NonNegativeReals,initialize=300)
m.TB_tb_ave = Var(m.t,within=NonNegativeReals,initialize=375)
m.TB_tb_in = Var(m.t,within=NonNegativeReals,initialize=460)
m.TB_tb_out = Var(m.t,within=NonNegativeReals,initialize=350)
m.dH_IPA = Var(within=NonNegativeReals)
m.dH_CO2 = Var(within=NonNegativeReals)
m.dH_B = Var(m.t,within=NonNegativeReals) # Latent heat of evap of liquid from stripper 
m.Pvap_B_CO2 = Var(m.t,within=NonNegativeReals)
m.Pvap_B_IPA = Var(m.t,within=NonNegativeReals,initialize=0.84)
m.P_B = Var(within=NonNegativeReals,initialize=40)
m.F_B = Var(m.t,within=NonNegativeReals,initialize=3.0)
m.F_B_bypass = Var(m.t,within=NonNegativeReals,initialize=4.5)

m.TB_tb_avedot = DerivativeVar(m.TB_tb_ave,wrt=m.t,initialize=0)

# Cooler:
m.rhoC_CO2 = Var(m.t,within=NonNegativeReals)
m.CpC_CO2 = Var(m.t,within=NonNegativeReals,initialize=40)
m.CpC_H2O = Var(m.t,within=NonNegativeReals,initialize=36.9)
# tube=solvent, shell=coolant
m.TC_tb_in = Var(m.t,within=NonNegativeReals,initialize=424)
m.TC_tb_ave = Var(m.t,within=NonNegativeReals,initialize=369)
m.TC_tb_out = Var(m.t,within=NonNegativeReals,initialize=313)
m.TC_sh_in = Var(m.t,within=NonNegativeReals,initialize=293)
m.TC_sh_out = Var(m.t,within=NonNegativeReals,initialize=306)
m.TC_sh_ave = Var(m.t,within=NonNegativeReals,initialize=299)
m.F_C = Var(m.t,within=NonNegativeReals,initialize=7.7)
m.F_C_H2O = Var(m.t,within=NonNegativeReals,initialize=22.25)
m.yC_IPA = Var(m.t,within=NonNegativeReals,initialize=0.0005)
m.yC_CO2 = Var(m.t,within=NonNegativeReals,initialize=1)

m.TC_tb_avedot = DerivativeVar(m.TC_tb_ave,wrt=m.t,initialize=0)
m.TC_sh_avedot = DerivativeVar(m.TC_sh_ave,wrt=m.t,initialize=0)

        
#########################
### Model Constraints ###
#########################

m.t.pprint()

# Extractor:
# inputs to this unit are solvent (FSE) and H2O/IPA "feed" (FSE)

# Differential:
def const_E1_rule(m,t):
    return m.ME*m.xE_IPAdot[m.NTE,t] == m.FFE[t]*(m.xE_F_IPA[t] - m.xE_IPA[m.NTE,t]) + \
               m.FSE[t]*(m.yE_IPA[m.NTE-1,t] - m.yE_IPA[m.NTE,t])
m.const_E1 = Constraint(m.t,rule=const_E1_rule)

def const_E2_rule(m,i,t):
    return m.ME*m.xE_IPAdot[i,t] == m.FFE[t]*(m.xE_IPA[i+1,t] - m.xE_IPA[i,t]) + \
               m.FSE[t]*(m.yE_IPA[i-1,t] - m.yE_IPA[i,t]) 
m.const_E2 = Constraint(m.set_E_int,m.t,rule=const_E2_rule)

def const_E3_rule(m,t):
    return m.ME*m.xE_IPAdot[1,t] == m.FFE[t]*(m.xE_IPA[2,t] - m.xE_IPA[1,t]) + \
               m.FSE[t]*(m.yE_S_IPA[t] - m.yE_IPA[1,t]) 
m.const_E3 = Constraint(m.t,rule=const_E3_rule)

# Algebraic
# vapor liquid equilium 
def const_E4_rule(m,i,t):
    return m.KE_IPA[t]*m.xE_IPA[i,t]*m.EMVE == m.yE_IPA[i,t]-m.yE_IPA[i-1,t]*(1-m.EMVE)
m.const_E4 = Constraint(m.set_E_2_NT,m.t,rule=const_E4_rule)

def const_E5_rule(m,t): 
    return m.KE_IPA[t]*m.xE_IPA[1,t]*m.EMVE == m.yE_IPA[1,t]-m.yE_S_IPA[t]*(1-m.EMVE)
m.const_E5 = Constraint(m.t,rule=const_E5_rule)

def const_E6_rule(m,t):
    return m.KE_IPA[t] == -0.002*m.TE[t] + 1.016
m.const_E6 = Constraint(m.t,rule=const_E6_rule)

# Feed comp disturbance constraint
def const_E7_rule(m,t):
    return m.xE_F_IPA[t] == m.W4_xF[t]
m.const_E7 = Constraint(m.t,rule=const_E7_rule)

# Acting as connectivity constraint 
def const_E8_rule(m,t):
    return m.TE[t] == m.TC_tb_out[t]
m.const_E8 = Constraint(m.t,rule=const_E8_rule)

# Feed flow disturbance constraint
def const_E9_rule(m,t):
    return m.FFE[t] == m.W3_F[t]
m.const_E9 = Constraint(m.t,rule=const_E9_rule)

# this constraint is deactivated because flow of solvent to extractor should equal 
# flow rate out of the cooler, always, and only necessarily must equal the nominal value of 
# 7.75 kmol/hr at t=0 (this initial condition still must be specified somewhere...)
#def const_E10_rule(m,t):
#    return m.FSE[t] == m.FSE_init
#m.const_E10 = Constraint(m.t,rule=const_E10_rule)

#def const_E11_rule(m):
#    return m.yE_S_IPA == m.yE_S_sp
#m.const_E11 = Constraint(rule=const_E11_rule)


# Flash:

# Equilibrium and split fraction
# why are there two VLE coefficients for flash
# and only one for the extractor? 

def const_F1_rule(m,t):
    return m.yF_CO2[t] == m.KF_CO2*m.xF_CO2[t]
m.const_F1 = Constraint(m.t,rule=const_F1_rule)

def const_F2_rule(m,t): 
    return (1-m.yF_CO2[t]) == m.KF_IPA*(1-m.xF_CO2[t])
m.const_F2 = Constraint(m.t,rule=const_F2_rule) 

def const_F3_rule(m,t): 
    return (1-m.yE_IPA[m.NTE.value,t]) == m.qF[t]*m.yF_CO2[t] + (1-m.qF[t])*m.xF_CO2[t] 
m.const_F3 = Constraint(m.t,rule=const_F3_rule) 


# Operating constraints

def const_F4_rule(m):
    return m.KF_CO2 == 1.3 #1.585
m.const_F4 = Constraint(rule=const_F4_rule)

def const_F5_rule(m):
    return m.KF_IPA == 0.0125
m.const_F5 = Constraint(rule=const_F5_rule) 
# These constants obtained from ... ?

#def const_F14_rule(m):
#    return m.TF == 313
#m.const_F14 = Constraint(rule=const_F14_rule)

# Thermo Equations

def const_F6_rule(m,t):
    A =  25.00
    B =  55.19
    C = -33.69
    D =  7.948
    E = -0.1366
    return m.Cp_CO2_v[t] == A+ B*m.TF[t]/1000 + C*(m.TF[t]/1000)**2 + D*(m.TF[t]/1000)**3 + E/(m.TF[t]/1000)**2
m.const_F6 = Constraint(m.t,rule=const_F6_rule)

def const_F7_rule(m,t):
    A =  25.00
    B =  55.19
    C = -33.69
    D =  7.948
    E = -0.1366
    F = -403.61
    H = -393.52
    return m.H_CO2_v[t] == A*(m.TF[t]/1000)+B/2*(m.TF[t]/1000)**2+C/3*(m.TF[t]/1000)**3+D/4*(m.TF[t]/1000)**4-E/(m.TF[t]/1000)+F-H
m.const_F7 = Constraint(m.t,rule=const_F7_rule)

def const_F8_rule(m,t):
    return m.Cp_CO2_l[t] == 195.67 + (m.TF[t]-310)*10.18
m.const_F8 = Constraint(m.t,rule=const_F8_rule)


def const_F9_rule(m,t): 
    H_CO2_l_298 = ((195.67-310*10.18)*298.15 + 10.18/2*298.15**2)/1000
    return m.H_CO2_l[t] == ((195.67-310*10.18)*m.TF[t] + 10.18/2*m.TF[t]**2)/1000 - H_CO2_l_298
m.const_F9 = Constraint(m.t,rule=const_F9_rule)

def const_F10_rule(m,t):
    return m.Cp_IPA_v[t] == 89.74+ (m.TF[t]-300)*0.219
m.const_F10 = Constraint(m.t,rule=const_F10_rule)

def const_F11_rule(m,t):
    H_IPA_v_298 = ((89.74-300*0.219)*298.15 + 0.219/2*298.15**2)/1000
    return m.H_IPA_v[t] == ((89.74-300*0.219)*m.TF[t] + 0.219/2*m.TF[t]**2)/1000 - H_IPA_v_298
m.const_F11 = Constraint(m.t,rule=const_F11_rule)

def const_F12_rule(m,t):
    return m.Cp_IPA_l[t] == 165.6+(m.TF[t]-311.6)*0.756
m.const_F12 = Constraint(m.t,rule=const_F12_rule)

def const_F13_rule(m,t):
    H_IPA_l_298 = (165.6-311.6*0.756)*298.15 + 0.756/2*298.15**2
    return m.H_IPA_l[t] == ((165.6-311.6*0.756)*m.TF[t] + 0.756/2*m.TF[t]**2-H_IPA_l_298)/1000
m.const_F13 = Constraint(m.t,rule=const_F13_rule)

# Pre-flash enthalpies of IPA and CO2 (using temp from extractor)

# treating sc phase as a liquid:
def const_F14_rule(m,t):
    H_CO2_l_298 = ((195.67-310*10.18)*298.15 + 10.18/2*298.15**2)/1000
    return m.H_CO2_E[t] == ((195.67-310*10.18)*m.TE[t] + 10.18/2*m.TE[t]**2)/1000 - H_CO2_l_298
m.const_F14 = Constraint(m.t,rule=const_F14_rule)

def const_F15_rule(m,t):
    H_IPA_l_298 = (165.6-311.6*0.756)*298.15 + 0.756/2*298.15**2
    return m.H_IPA_E[t] == ((165.6-311.6*0.756)*m.TE[t] + 0.756/2*m.TE[t]**2-H_IPA_l_298)/1000
m.const_F15 = Constraint(m.t,rule=const_F15_rule)

# treating sc phase as a vapor:
#def const_F14_rule(m):
#    A =  25.00
#    B =  55.19
#    C = -33.69
#    D =  7.948
#    E = -0.1366
#    F = -403.61
#    H = -393.52
#    return m.H_CO2_E == A*(m.TE/1000)+B/2*(m.TE/1000)**2+C/3*(m.TE/1000)**3+D/4*(m.TE/1000)**4-E/(m.TE/1000)+F-H
#m.const_F14 = Constraint(rule=const_F14_rule)
#
#def const_F15_rule(m): 
#    H_IPA_v_298 = ((89.74-300*0.219)*298.15 + 0.219/2*298.15**2)/1000
#    return m.H_IPA_E == ((89.74-300*0.219)*m.TE + 0.219/2*m.TE**2)/1000 - H_IPA_v_298
#m.const_F15 = Constraint(rule=const_F15_rule)

# Enthalpy balance, used to calculate T:

def const_F16_rule(m,t):
    return m.yE_IPA[m.NTE.value,t]*m.H_IPA_E[t] + (1-m.yE_IPA[m.NTE.value,t])*m.H_CO2_E[t] == \
            (1-m.qF[t])*(m.xF_CO2[t]*m.H_CO2_l[t] + (1-m.xF_CO2[t])*m.H_IPA_l[t]) + \
            m.qF[t]*(m.yF_CO2[t]*m.H_CO2_v[t] + (1-m.yF_CO2[t])*m.H_IPA_v[t])
m.const_F16 = Constraint(m.t,rule=const_F16_rule)

# Next: use the outlet from the flash (x/yF_CO2) to specify inlet to the Stripper
# Q: should I specify x/yF_IPA as well, or just use CO2...
# probably should specify

def const_F17_rule(m,t):
    return m.yF_IPA[t] == 1-m.yF_CO2[t]
m.const_F17 = Constraint(m.t,rule=const_F17_rule)

def const_F18_rule(m,t):
    return m.xF_IPA[t] == 1-m.xF_CO2[t]
m.const_F18 = Constraint(m.t,rule=const_F18_rule)

# connectivity constraint, extractor to flash
def const_F19_rule(m,t):
    return m.FFF[t] == m.FSE[t]
m.const_F19 = Constraint(m.t,rule=const_F19_rule)

def const_F20_rule(m):
    return m.P1_F == 100 # atm
m.const_F20 = Constraint(rule=const_F20_rule)

def const_F21_rule(m):
    return m.P2_F == 40 # atm
m.const_F22 = Constraint(rule=const_F21_rule)

# Stripping column
# will need to write my own constraints to 
# keep my variables organized... 

# operating constraints:
# (defining nominal values of manipulated variables)

def const_S1_rule(m,t):
    return m.FMK[t] == m.W1_MK[t]
m.const_S1 = Constraint(m.t,rule=const_S1_rule)

### control inputs * : should change these constraints (add U1/U2 parameters w/ 1 col pt...)
def const_S2_rule(m,t):
    return m.FV_S[t] == m.U1_V[t]
m.const_S2 = Constraint(m.t,rule=const_S2_rule)

def const_S3_rule(m,t):
    return m.FL_S[t] == m.U2_L[t]
m.const_S3 = Constraint(m.t,rule=const_S3_rule)

#def const_S4_rule(m):
#    return m.FSS == m.FSS_nom
#m.const_S4 = Constraint(rule=const_S4_rule)

# connectivity constraints:
def const_S5_rule(m,t):
    return m.FF_S[t] == m.FFF[t]
m.const_S5 = Constraint(m.t,rule=const_S5_rule)

def const_S6_rule(m,t):
    return m.xF_S[t] == m.qF[t]*m.yF_CO2[t] + (1-m.qF[t])*m.xF_CO2[t]
m.const_S6 = Constraint(m.t,rule=const_S6_rule)

# differential
# (mass balance)

# No "differential" mass balance for top tray (condenser)

# tray 4:
def const_S7_rule(m,t): 
    # need to re-calculate xS_CO2[NTS]
    return m.MS*m.xS_CO2dot[m.NTS-1,t] == \
                                  m.FL_S[t]*(m.xS_CO2[m.NTS.value,t] - m.xS_CO2[4,t]) + \
            (m.FV_S[t] + m.qF[t]*m.FF_S[t])*(m.yS_CO2[3,t] - m.yS_CO2[4,t])
m.const_S7 = Constraint(m.t,rule=const_S7_rule)

# tray 3 (feed):
def const_S8_rule(m,t):
    return m.MS*m.xS_CO2dot[m.NTFS,t] == \
                                  m.FL_S[t]*(m.xS_CO2[4,t] - m.xS_CO2[m.NTFS.value,t]) + \
                                  m.FV_S[t]*(m.yS_CO2[2,t] - m.xS_CO2[m.NTFS.value,t]) + \
            (m.FV_S[t] + m.qF[t]*m.FF_S[t])*(m.xS_CO2[m.NTFS.value,t] - m.yS_CO2[m.NTFS.value,t]) + \
                                  m.FF_S[t]*(m.xF_S[t] - m.xS_CO2[m.NTFS.value,t]) 
m.const_S8 = Constraint(m.t,rule=const_S8_rule)

# tray 2:
def const_S9_rule(m,t):
    return m.MS*m.xS_CO2dot[2,t] == \
            (m.FL_S[t] + (1-m.qF[t])*m.FF_S[t])*(m.xS_CO2[3,t] - m.xS_CO2[2,t]) + \
                                      m.FV_S[t]*(m.yS_CO2[1,t] - m.yS_CO2[2,t])
m.const_S9 = Constraint(m.t,rule=const_S9_rule)

# tray 1 (reboiler):
def const_S10_rule(m,t):
    return m.MS*m.xS_CO2dot[1,t] == \
            (m.FL_S[t] + (1-m.qF[t])*m.FF_S[t])*(m.xS_CO2[2,t] - m.xS_CO2[1,t]) + \
                                      m.FV_S[t]*(m.xS_CO2[1,t] - m.yS_CO2[1,t])
m.const_S10 = Constraint(m.t,rule=const_S10_rule)

# algebraic
# vapor liquid (thermal...) equilibrium
# perfect control in condenser/boiler

# condenser equations:
#def const_S11_rule(m):
#    return (m.FVS+m.qF*m.FFS)*m.yS_CO2[4] + m.FMK == (m.FFS+m.FLS)*m.yS_CO2[5]
#m.const_S11 = Constraint(rule=const_S11_rule)

def const_S19_rule(m,t):
    return m.FT_S[t] == m.FV_S[t] + m.qF[t]*m.FF_S[t]
m.const_S19 = Constraint(m.t,rule=const_S19_rule)

def const_S20_rule(m,t):
    return m.F_comp[t] == m.FMK[t] + m.FT_S[t]
m.const_S20 = Constraint(m.t,rule=const_S20_rule)

def const_S11_rule(m,t):
    # This equation specifies composition in "condenser," after 
    # mixing with make-up CO2 stream
    return m.FT_S[t]*m.yS_CO2[4,t] + m.FMK[t] == m.F_comp[t]*m.yS_CO2[5,t]
m.const_S11 = Constraint(m.t,rule=const_S11_rule)

#def const_S12_rule(m): 
    # This equation specifies FSS var
#    return m.FVS+m.qF*m.FFS+m.FMK == m.FLS+m.FSS
#m.const_S12 = Constraint(rule=const_S12_rule)

def const_S12_rule(m,t):
    # Flow in the compressor equals combined reflux and distillate flow rates
    return m.F_comp[t] == m.FL_S[t] + m.FD_S[t]
m.const_S12 = Constraint(m.t,rule=const_S12_rule)

# reboiler bottoms equation:
def const_S13_rule(m,t):
    return m.FB_S[t] == m.FL_S[t] + (1-m.qF[t])*m.FF_S[t] - m.FV_S[t]
m.const_S13 = Constraint(m.t,rule=const_S13_rule)

# VLE:
def const_S14_rule(m):
    return m.KS_CO2 == 1.3
m.const_S14 = Constraint(rule=const_S14_rule)

def const_S18_rule(m):
    return m.KS_IPA == 0.0213
m.const_S18 = Constraint(rule=const_S18_rule)

def const_S15_rule(m,i,t):
    return m.KS_IPA*m.EMVS*(1-m.xS_CO2[i,t]) == \
            (1-m.yS_CO2[i,t]) - (1-m.yS_CO2[i-1,t])*(1-m.EMVS)
m.const_S15 = Constraint(m.set_S_eq,m.t,rule=const_S15_rule)

# assume reboiler reaches full vapor liquid equilibrium:
def const_S16_rule(m,t):
    return m.KS_IPA*(1-m.xS_CO2[1,t]) == (1-m.yS_CO2[1,t])
m.const_S16 = Constraint(m.t,rule=const_S16_rule)

# compositions are the same across phases in the condenser
# (because they are different 'branches' of the same stream)
def const_S17_rule(m,t): 
    return m.yS_CO2[m.NTS.value,t] == m.xS_CO2[m.NTS.value,t]
m.const_S17 = Constraint(m.t,rule=const_S17_rule)

# Condenser:
def const_cond1_rule(m):
    return m.P1_cond == m.P2_F
m.const_cond1 = Constraint(rule=const_cond1_rule)
    
def const_cond2_rule(m,t):
    return m.T1_cond[t] == m.TF[t]
m.const_cond2 = Constraint(m.t,rule=const_cond2_rule)

def const_cond3_rule(m,t):
    return m.V1_cond[t] == (1/m.P1_cond)*m.F_comp[t]*0.08206*m.T1_cond[t]
m.const_cond3 = Constraint(m.t,rule=const_cond3_rule)

def const_cond4_rule(m):
    return m.P2_cond == 100.0
m.const_cond4 = Constraint(rule=const_cond4_rule)

def const_cond5_rule(m,t):
    return m.P1_cond*m.V1_cond[t]**1.4 == m.P2_cond*m.V2_cond[t]**1.4
m.const_cond5 = Constraint(m.t,rule=const_cond5_rule)

def const_cond6_rule(m,t):
    return m.P2_cond*m.V2_cond[t] == m.F_comp[t]*0.08206*m.T2_cond[t]
m.const_cond6 = Constraint(m.t,rule=const_cond6_rule)


# Reboiler
# reboiler already almost specified, just need energy balance
# need to:
#   enforce (mixture) heat capacity as a function of temp (actually not necessary for st.st.)
#   maybe enforce latent heat (of vaporization) of liquid from stripper (x2/L2?) (as fcn of Temp?)
#   calculate reboiler temperature
#   calculate temperature of vapor (solvent) leaving reboiler (with steadied en. bal. constraint)
#   ? calculate density of vapor CO2 solvent through reboiler ? (Ramchandran uses average)
#   ^ again not necessary for steady state

def const_B1_rule(m,t):
    return m.Pvap_B_IPA[t]*101.325 == 10**(6.866 - 1360/(197.6+m.TB_sh[t]-273.15))
m.const_B1 = Constraint(m.t,rule=const_B1_rule)

def const_B2_rule(m):
    return m.P_B == 40.0
m.const_B2 = Constraint(rule=const_B2_rule)

def const_B3_rule(m,t):
    return m.Pvap_B_IPA[t]*(1-m.xS_CO2[1,t]) == m.P_B*(1-m.yS_CO2[1,t])
m.const_B3 = Constraint(m.t,rule=const_B3_rule)

def const_B4_rule(m):
    return m.dH_IPA == 44000 # kJ/kmol
m.const_B4 = Constraint(rule=const_B4_rule)

def const_B5_rule(m):
    return m.dH_CO2 == 15326 # kJ/kmol
m.const_B5 = Constraint(rule=const_B5_rule)

def const_B6_rule(m,t):
    return m.dH_B[t] == m.yS_CO2[1,t]*m.dH_CO2 + (1-m.yS_CO2[1,t])*m.dH_IPA
m.const_B6 = Constraint(m.t,rule=const_B6_rule)

# energy balance for temperature dynamics, really should be a PDE...
def const_B7_rule(m,t):
    return m.F_B[t]*m.CpB_CO2[t]*m.TB_tb_avedot[t] == \
            m.CpB_CO2[t]*m.F_B[t]*(m.TB_tb_in[t]-m.TB_tb_out[t]) \
            + m.U_ov_B*m.A_B*(m.TB_sh[t] - m.TB_tb_ave[t])
m.const_B7 = Constraint(m.t,rule=const_B7_rule)

# constant heat capacity, that of CO2, for simplicity
def const_B8_rule(m,t):
    return m.CpB_CO2[t] == 40.0
m.const_B8 = Constraint(m.t,rule=const_B8_rule)

def const_B9_rule(m,t):
    return m.TB_tb_in[t] == m.T2_cond[t]
m.const_B9 = Constraint(m.t,rule=const_B9_rule)

# is this constraint correct, or should VdH=UAdT ? 
def const_B10_rule(m,t):
    return (m.TB_tb_in[t] - m.TB_tb_out[t])*m.CpB_CO2[t]*m.F_B[t] == m.FV_S[t]*m.dH_B[t]
m.const_B10 = Constraint(m.t,rule=const_B10_rule)

def const_B11_rule(m,t):
    return 2*m.TB_tb_ave[t] - m.TB_tb_in[t] == m.TB_tb_out[t]
m.const_B11 = Constraint(m.t,rule=const_B11_rule)

def const_B12_rule(m,t):
    return m.F_B_bypass[t] == m.FD_S[t] - m.F_B[t]
m.const_B12 = Constraint(m.t,rule=const_B12_rule)

# Cooler

def const_C1_rule(m,t):
    return 2*m.TC_tb_ave[t] == m.TC_tb_out[t] + m.TC_tb_in[t]
m.const_C1 = Constraint(m.t,rule=const_C1_rule)

def const_C2_rule(m,t):
    return m.F_C[t] == m.F_B[t] + m.F_B_bypass[t]
m.const_C2 = Constraint(m.t,rule=const_C2_rule)

def const_C3_rule(m,t):
    return m.TC_tb_in[t]*m.F_C[t] == m.TB_tb_out[t]*m.F_B[t] + m.T2_cond[t]*m.F_B_bypass[t]
m.const_C3 = Constraint(m.t,rule=const_C3_rule)

#def const_C4_rule(m):
#    return m.TC_cool == 306 
#m.const_C4 = Constraint(rule=const_C4_rule)

def const_C5_rule(m,t):
    return m.CpC_CO2[t] == m.CpB_CO2[t]
m.const_C5 = Constraint(m.t,rule=const_C5_rule)

# assuming heat capacity of stream is that of CO2
# Really should be a PDE...
def const_C6_rule(m,t):
    return m.CpC_CO2[t]*m.F_C[t]*m.TC_tb_avedot[t] == \
            m.F_C[t]*m.CpC_CO2[t]*(m.TC_tb_in[t] - m.TC_tb_out[t]) + \
            m.U_ov_C*m.A_C*       (m.TC_sh_ave[t] - m.TC_tb_ave[t])
m.const_C6 = Constraint(m.t,rule=const_C6_rule)

def const_C7_rule(m,t):
    return m.CpC_H2O[t] == 36.9
m.const_C7 = Constraint(m.t,rule=const_C7_rule)

def const_C8_rule(m,t):
    return m.TC_sh_in[t] == m.W2_Tcool[t]
m.const_C8 = Constraint(m.t,rule=const_C8_rule)

def const_C9_rule(m,t):
    return 2*m.TC_sh_ave[t] == m.TC_sh_in[t] + m.TC_sh_out[t]
m.const_C9 = Constraint(m.t,rule=const_C9_rule)

def const_C10_rule(m,t):
    return m.F_C_H2O[t] == m.U3_Fcool[t]
m.const_C10 = Constraint(m.t,rule=const_C10_rule)

# Really should be a PDE...
def const_C11_rule(m,t):
    return m.CpC_H2O[t]*m.F_C_H2O[t]*m.TC_sh_avedot[t] == \
            m.F_C_H2O[t]*m.CpC_H2O[t]*(m.TC_sh_in[t] - m.TC_sh_out[t]) + \
            m.U_ov_C*m.A_C           *(m.TC_tb_ave[t] - m.TC_sh_ave[t])
m.const_C11 = Constraint(m.t,rule=const_C11_rule)

def const_C12_rule(m,t):
    return m.yC_CO2[t] == m.xS_CO2[m.NTS,t]
m.const_C12 = Constraint(m.t,rule=const_C12_rule)

def const_C13_rule(m,t):
    return m.yC_IPA[t] == 1 - m.yC_CO2[t]
m.const_C13 = Constraint(m.t,rule=const_C13_rule)

# connectivity
def const_C14_rule(m,t):
    return m.yC_IPA[t] == m.yE_S_IPA[t]
m.const_C14 = Constraint(m.t,rule=const_C14_rule)

#...not a chance this converges...
def const_C15_rule(m,t):
    return m.F_C[t] == m.FSE[t]
    #return m.FSE == m.FMK + m.FT_S - m.FD_S
    #return m.FMK + m.FV_S == m.FL_S
    #return m.FMK == m.FL_S-m.FV_S+(1-m.qF)*m.FF_S
m.const_C15 = Constraint(m.t,rule=const_C15_rule)

#####################################################
### Initial Conditions for Differential Variables ###
#####################################################

def const_init_xE_rule(m,i):
    return m.xE_IPA[i,0] == m.xE_IPA_0[i]
m.const_init_xE = Constraint(m.set_E_Tray,rule=const_init_xE_rule)

def const_init_xS_rule(m,i):
    return m.xS_CO2[i,0] == m.xS_CO2_0[i]
m.const_init_xS = Constraint(m.set_S_dyn,rule=const_init_xS_rule)

def const_init_TB_tb_rule(m):
    return m.TB_tb_ave[0] == m.TB_tb_0
m.const_init_TB_tb = Constraint(rule=const_init_TB_tb_rule)

def const_init_TC_tb_rule(m):
    return m.TC_tb_ave[0] == m.TC_tb_0
m.const_init_TC_tb = Constraint(rule=const_init_TC_tb_rule)

def const_init_TC_sh_rule(m):
    return m.TC_sh_ave[0] == m.TC_sh_0
m.const_init_TC_sh = Constraint(rule=const_init_TC_sh_rule)

########################
### Discretize Model ###
########################

discretizer = TransformationFactory('dae.collocation')
discretizer.apply_to(m,nfe=25,ncp=4,scheme='LAGRANGE-RADAU')

# Reduce number of collocation points for control variables 
# to force them to be piecewise constant 

#m = discretizer.reduce_collocation_points(m,var=m.U1_V,ncp=1,contset=m.t)
#m = discretizer.reduce_collocation_points(m,var=m.U2_L,ncp=1,contset=m.t)
#m = discretizer.reduce_collocation_points(m,var=m.U3_Fcool,ncp=1,contset=m.t)


###################
### objective #####
###################

#m.obj_ss = Objective(sense=minimize,expr=(m.TC_tb_out-313.0)**2)
# needs to be adjusted to sum over all time? Or not because I have degrees of freedom now...?
#m.obj_ss = Objective(sense=minimize,expr=(m.F_C-m.FSE)**2+(m.TC_tb_out-313.0)**2)
m.obj_ss = Objective(sense=minimize,expr=0.0)
