from pyomo.environ import * 
from pyomo.dae import ContinuousSet, DerivativeVar
import numpy as np

m = ConcreteModel()

# Control options

# Model options
const_Cp_switch = 0

# Discretization options

m.t = ContinuousSet(bounds=(0,300)) # solve nmpc opt problem over 5 hours

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
print('A_B: \t'+str(m.A_B.value))
print('A_cr: \t' + str(m.A_cr_B.value))
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
m.set_S_eq.pprint()

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

#######################
### State variables ###
#######################

# Extractor:
# compositions below are that of IPA 
m.xE_IPA = Var(m.set_E_Tray,within=NonNegativeReals,initialize=0.1)
m.yE_IPA = Var(m.set_E_Tray,within=NonNegativeReals,initialize=0.1)
m.xE_F_IPA = Var(within=NonNegativeReals,initialize=0.19)
m.yE_S_IPA = Var(within=NonNegativeReals,initialize=0.0002)
m.FFE = Var(within=NonNegativeReals,initialize=1.8) # Feed Flow into Extractor
# ^is this a parameter, state variable, or control variable? 
m.FSE = Var(within=NonNegativeReals,initialize=7.75) # Solvent Flow into Extractor
m.KE_IPA = Var(within=NonNegativeReals,initialize=1.3)
m.TE = Var(within=NonNegativeReals,initialize=313)


# Flash:
m.FFF = Var(within=NonNegativeReals,initialize=m.FSE_init) # Flow of Feed to Flash
m.qF = Var(within=NonNegativeReals,initialize=0.999) # Vapor split fraction after flash
m.xF_IPA = Var(within=NonNegativeReals,initialize=0.2)
m.xF_CO2 = Var(within=NonNegativeReals,initialize=0.9)
m.yF_IPA = Var(within=NonNegativeReals,initialize=0.004)
m.yF_CO2 = Var(within=NonNegativeReals,initialize=0.9)
m.KF_CO2 = Var(within=NonNegativeReals,initialize=1.5)
m.KF_IPA = Var(within=NonNegativeReals,initialize=0.02)
m.Cp_CO2_l = Var(within=NonNegativeReals,initialize=226.21)
m.Cp_CO2_v = Var(within=NonNegativeReals,initialize=37.82)
m.Cp_IPA_l = Var(within=NonNegativeReals,initialize=166.66)
m.Cp_IPA_v = Var(within=NonNegativeReals,initialize=92.59)
m.H_CO2_l = Var(within=Reals,initialize=3.22) # why are the liquid streams at higher enthalpy?
m.H_CO2_v = Var(within=Reals,initialize=0.55)
m.H_CO2_E = Var(within=Reals,initialize=3.0)
m.H_IPA_l = Var(within=Reals,initialize=2.383)
m.H_IPA_v = Var(within=Reals,initialize=1.34)
m.H_IPA_E = Var(within=Reals,initialize=2.0)
m.TF = Var(within=NonNegativeReals,initialize=310)
m.P1_F = Var(within=NonNegativeReals,initialize=100)
m.P2_F = Var(within=NonNegativeReals,initialize=40)

# Stripping Column:
m.xS_CO2 = Var(m.set_S_Tray,within=NonNegativeReals,initialize=0.9)
m.yS_CO2 = Var(m.set_S_Tray,within=NonNegativeReals,initialize=0.9)
m.FF_S = Var(within=NonNegativeReals,initialize=m.FSE_init) # feed flow rate
m.FL_S = Var(within=NonNegativeReals,initialize=m.FLS_nom.value) # liquid reflux flow rate
m.FV_S = Var(within=NonNegativeReals,initialize=m.FVS_nom.value) # vapor boilup flow rate
m.FB_S = Var(within=NonNegativeReals,initialize=m.FMK_nom.value) # bottoms product flow rate 
m.FMK = Var(within=NonNegativeReals,initialize=m.FMK_nom.value) # make up CO2 flow rate
#m.FS_S = Var(within=NonNegativeReals,initialize=m.FSS_nom.value) 
m.FT_S = Var(within=NonNegativeReals,initialize=m.FTS_nom.value) # flow off top of column
m.FD_S = Var(within=NonNegativeReals,initialize=m.FDS_nom.value) # distillate (solvent) flow rate
m.F_comp = Var(within=NonNegativeReals,initialize=m.F_comp_nom.value) # flow rate through compresser 
m.xF_S = Var(within=NonNegativeReals,initialize=0.9) # CO2 fraction of combined (liq & vap) feed to stripper
m.KS_CO2 = Var(within=NonNegativeReals,initialize=1.3)
m.KS_IPA = Var(within=NonNegativeReals,initialize=m.KS_IPA_nom.value)
m.TS = Var(within=NonNegativeReals,initialize=350)
m.PS = Var(within=NonNegativeReals,initialize=40)

# Condenser:
m.P1_cond = Var(within=NonNegativeReals,initialize=40)
m.P2_cond = Var(within=NonNegativeReals,initialize=100)
m.T1_cond = Var(within=NonNegativeReals,initialize=350)
m.T2_cond = Var(within=NonNegativeReals,initialize=460)
m.V1_cond = Var(within=NonNegativeReals,initialize=8)
m.V2_cond = Var(within=NonNegativeReals,initialize=3)

# Reboiler:
#m.rhoB_CO2 = Var(within=NonNegativeReals) # Density in kg/m^3
m.CpB_CO2 = Var(within=NonNegativeReals,initialize=40)
m.TB_sh = Var(within=NonNegativeReals,initialize=300)
m.TB_tb_ave = Var(within=NonNegativeReals,initialize=375)
m.TB_tb_in = Var(within=NonNegativeReals,initialize=460)
m.TB_tb_out = Var(within=NonNegativeReals,initialize=350)
m.dH_IPA = Var(within=NonNegativeReals,initialize=44000)
m.dH_CO2 = Var(within=NonNegativeReals,initialize=15326)
m.dH_B = Var(within=NonNegativeReals,initialize=15500) # Latent heat of evap of liquid from stripper 
#m.Pvap_B_CO2 = Var(within=NonNegativeReals)
m.Pvap_B_IPA = Var(within=NonNegativeReals,initialize=0.84)
m.P_B = Var(within=NonNegativeReals,initialize=40)
m.F_B = Var(within=NonNegativeReals,initialize=3.0)
m.F_B_bypass = Var(within=NonNegativeReals,initialize=4.5)

# Cooler:
#m.rhoC_CO2 = Var(within=NonNegativeReals)
m.CpC_CO2 = Var(within=NonNegativeReals,initialize=40)
m.CpC_H2O = Var(within=NonNegativeReals,initialize=36.9)
# tube=solvent, shell=coolant
m.TC_tb_in = Var(within=NonNegativeReals,initialize=424)
m.TC_tb_ave = Var(within=NonNegativeReals,initialize=369)
m.TC_tb_out = Var(within=NonNegativeReals,initialize=313)
m.TC_sh_in = Var(within=NonNegativeReals,initialize=293)
m.TC_sh_out = Var(within=NonNegativeReals,initialize=306)
m.TC_sh_ave = Var(within=NonNegativeReals,initialize=299)
m.F_C = Var(within=NonNegativeReals,initialize=7.7)
m.F_C_H2O = Var(within=NonNegativeReals,initialize=22.25)
m.yC_IPA = Var(within=NonNegativeReals,initialize=0.0005)
m.yC_CO2 = Var(within=NonNegativeReals,initialize=1)

        
#########################
### Model Constraints ###
#########################

# Extractor:
# inputs to this unit are solvent (FSE) and H2O/IPA "feed" (FSE)

# Differential:
def const_E1_rule(m):
    return 0 == m.FFE*(m.xE_F_IPA - m.xE_IPA[m.NTE]) + \
               m.FSE*(m.yE_IPA[m.NTE-1] - m.yE_IPA[m.NTE])
m.const_E1 = Constraint(rule=const_E1_rule)

def const_E2_rule(m,i):
    return 0 == m.FFE*(m.xE_IPA[i+1] - m.xE_IPA[i]) + \
               m.FSE*(m.yE_IPA[i-1] - m.yE_IPA[i]) 
m.const_E2 = Constraint(m.set_E_int,rule=const_E2_rule)

def const_E3_rule(m):
    return 0 == m.FFE*(m.xE_IPA[2] - m.xE_IPA[1]) + \
               m.FSE*(m.yE_S_IPA - m.yE_IPA[1]) 
m.const_E3 = Constraint(rule=const_E3_rule)

# Algebraic
# vapor liquid equilium 
def const_E4_rule(m,i):
    return m.KE_IPA*m.xE_IPA[i]*m.EMVE == m.yE_IPA[i]-m.yE_IPA[i-1]*(1-m.EMVE)
m.const_E4 = Constraint(m.set_E_2_NT,rule=const_E4_rule)

def const_E5_rule(m): 
    return m.KE_IPA*m.xE_IPA[1]*m.EMVE == m.yE_IPA[1]-m.yE_S_IPA*(1-m.EMVE)
m.const_E5 = Constraint(rule=const_E5_rule)

def const_E6_rule(m):
    return m.KE_IPA == -0.002*m.TE + 1.016
m.const_E6 = Constraint(rule=const_E6_rule)

# Operating (initial) condition constraints
def const_E7_rule(m):
    return m.xE_F_IPA == m.xE_F_nom
m.const_E7 = Constraint(rule=const_E7_rule)

def const_E8_rule(m):
    return m.TE == m.TC_tb_out
    #return m.TE == m.TE_sp
m.const_E8 = Constraint(rule=const_E8_rule)

def const_E9_rule(m):
    return m.FFE == m.FFE_nom
m.const_E9 = Constraint(rule=const_E9_rule)

def const_E10_rule(m):
    return m.FSE == m.FSE_init
m.const_E10 = Constraint(rule=const_E10_rule)

#def const_E11_rule(m):
#    return m.yE_S_IPA == m.yE_S_sp
#m.const_E11 = Constraint(rule=const_E11_rule)


# Flash:

# Equilibrium and split fraction
# why are there two VLE coefficients for flash
# and only one for the extractor? 

def const_F1_rule(m):
    return m.yF_CO2 == m.KF_CO2*m.xF_CO2
m.const_F1 = Constraint(rule=const_F1_rule)

def const_F2_rule(m): 
    return (1-m.yF_CO2) == m.KF_IPA*(1-m.xF_CO2)
m.const_F2 = Constraint(rule=const_F2_rule) 

def const_F3_rule(m): 
    return (1-m.yE_IPA[m.NTE.value]) == m.qF*m.yF_CO2 + (1-m.qF)*m.xF_CO2 
m.const_F3 = Constraint(rule=const_F3_rule) 


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

def const_F6_rule(m):
    A =  25.00
    B =  55.19
    C = -33.69
    D =  7.948
    E = -0.1366
    return m.Cp_CO2_v == A+ B*m.TF/1000 + C*(m.TF/1000)**2 + D*(m.TF/1000)**3 + E/(m.TF/1000)**2
m.const_F6 = Constraint(rule=const_F6_rule)

def const_F7_rule(m):
    A =  25.00
    B =  55.19
    C = -33.69
    D =  7.948
    E = -0.1366
    F = -403.61
    H = -393.52
    return m.H_CO2_v == A*(m.TF/1000)+B/2*(m.TF/1000)**2+C/3*(m.TF/1000)**3+D/4*(m.TF/1000)**4-E/(m.TF/1000)+F-H
m.const_F7 = Constraint(rule=const_F7_rule)

def const_F8_rule(m):
    return m.Cp_CO2_l == 195.67 + (m.TF-310)*10.18
m.const_F8 = Constraint(rule=const_F8_rule)


def const_F9_rule(m): 
    H_CO2_l_298 = ((195.67-310*10.18)*298.15 + 10.18/2*298.15**2)/1000
    return m.H_CO2_l == ((195.67-310*10.18)*m.TF + 10.18/2*m.TF**2)/1000 - H_CO2_l_298
m.const_F9 = Constraint(rule=const_F9_rule)

def const_F10_rule(m):
    return m.Cp_IPA_v == 89.74+ (m.TF-300)*0.219
m.const_F10 = Constraint(rule=const_F10_rule)

def const_F11_rule(m):
    H_IPA_v_298 = ((89.74-300*0.219)*298.15 + 0.219/2*298.15**2)/1000
    return m.H_IPA_v == ((89.74-300*0.219)*m.TF + 0.219/2*m.TF**2)/1000 - H_IPA_v_298
m.const_F11 = Constraint(rule=const_F11_rule)

def const_F12_rule(m):
    return m.Cp_IPA_l == 165.6+(m.TF-311.6)*0.756
m.const_F12 = Constraint(rule=const_F12_rule)

def const_F13_rule(m):
    H_IPA_l_298 = (165.6-311.6*0.756)*298.15 + 0.756/2*298.15**2
    return m.H_IPA_l == ((165.6-311.6*0.756)*m.TF + 0.756/2*m.TF**2-H_IPA_l_298)/1000
m.const_F13 = Constraint(rule=const_F13_rule)

# Pre-flash enthalpies of IPA and CO2 (using temp from extractor)

# treating sc phase as a liquid:
def const_F14_rule(m):
    H_CO2_l_298 = ((195.67-310*10.18)*298.15 + 10.18/2*298.15**2)/1000
    return m.H_CO2_E == ((195.67-310*10.18)*m.TE + 10.18/2*m.TE**2)/1000 - H_CO2_l_298
m.const_F14 = Constraint(rule=const_F14_rule)

def const_F15_rule(m):
    H_IPA_l_298 = (165.6-311.6*0.756)*298.15 + 0.756/2*298.15**2
    return m.H_IPA_E == ((165.6-311.6*0.756)*m.TE + 0.756/2*m.TE**2-H_IPA_l_298)/1000
m.const_F15 = Constraint(rule=const_F15_rule)

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

def const_F16_rule(m):
    return m.yE_IPA[m.NTE.value]*m.H_IPA_E + (1-m.yE_IPA[m.NTE.value])*m.H_CO2_E == \
            (1-m.qF)*(m.xF_CO2*m.H_CO2_l + (1-m.xF_CO2)*m.H_IPA_l) + \
            m.qF*(m.yF_CO2*m.H_CO2_v + (1-m.yF_CO2)*m.H_IPA_v)
m.const_F16 = Constraint(rule=const_F16_rule)

# Next: use the outlet from the flash (x/yF_CO2) to specify inlet to the Stripper
# Q: should I specify x/yF_IPA as well, or just use CO2...
# probably should specify

def const_F17_rule(m):
    return m.yF_IPA == 1-m.yF_CO2
m.const_F17 = Constraint(rule=const_F17_rule)

def const_F18_rule(m):
    return m.xF_IPA == 1-m.xF_CO2
m.const_F18 = Constraint(rule=const_F18_rule)

# connectivity constraint, extractor to flash
def const_F19_rule(m):
    return m.FFF == m.FSE
m.const_F19 = Constraint(rule=const_F19_rule)

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
# (defining nominal calues of manipulated variables)
#def const_S1_rule(m):
#    return m.FMK == m.FMK_nom
#m.const_S1 = Constraint(rule=const_S1_rule)

def const_S2_rule(m):
    return m.FV_S == m.FVS_nom
m.const_S2 = Constraint(rule=const_S2_rule)

def const_S3_rule(m):
    return m.FL_S == m.FLS_nom
m.const_S3 = Constraint(rule=const_S3_rule)

#def const_S4_rule(m):
#    return m.FSS == m.FSS_nom
#m.const_S4 = Constraint(rule=const_S4_rule)

# connectivity constraints:
def const_S5_rule(m):
    return m.FF_S == m.FFF
m.const_S5 = Constraint(rule=const_S5_rule)

def const_S6_rule(m):
    return m.xF_S == m.qF*m.yF_CO2 + (1-m.qF)*m.xF_CO2
m.const_S6 = Constraint(rule=const_S6_rule)

# differential
# (mass balance)

# No "differential" mass balance for top tray (condenser)

# tray 4:
def const_S7_rule(m): 
    # need to re-calculate xS_CO2[NTS]
    return 0 == m.FL_S*(m.xS_CO2[m.NTS.value] - m.xS_CO2[4]) + \
            (m.FV_S + m.qF*m.FF_S)*(m.yS_CO2[3] - m.yS_CO2[4])
m.const_S7 = Constraint(rule=const_S7_rule)

# tray 3 (feed):
def const_S8_rule(m):
    return 0 == m.FL_S*(m.xS_CO2[4] - m.xS_CO2[m.NTFS.value]) + \
            m.FV_S*(m.yS_CO2[2] - m.xS_CO2[m.NTFS.value]) + \
            (m.FV_S + m.qF*m.FF_S)*(m.xS_CO2[m.NTFS.value] - m.yS_CO2[m.NTFS.value]) + \
            m.FF_S*(m.xF_S - m.xS_CO2[m.NTFS.value]) 
m.const_S8 = Constraint(rule=const_S8_rule)

# tray 2:
def const_S9_rule(m):
    return 0 == (m.FL_S+ (1-m.qF)*m.FF_S)*(m.xS_CO2[3] - m.xS_CO2[2]) + \
            m.FV_S*(m.yS_CO2[1] - m.yS_CO2[2])
m.const_S9 = Constraint(rule=const_S9_rule)

# tray 1 (reboiler):
def const_S10_rule(m):
    return 0 == (m.FL_S + (1-m.qF)*m.FF_S)*(m.xS_CO2[2]-m.xS_CO2[1]) + \
            m.FV_S*(m.xS_CO2[1] - m.yS_CO2[1])
m.const_S10 = Constraint(rule=const_S10_rule)

# algebraic
# vapor liquid (thermal...) equilibrium
# perfect control in condenser/boiler

# condenser equations:
#def const_S11_rule(m):
#    return (m.FVS+m.qF*m.FFS)*m.yS_CO2[4] + m.FMK == (m.FFS+m.FLS)*m.yS_CO2[5]
#m.const_S11 = Constraint(rule=const_S11_rule)

def const_S19_rule(m):
    return m.FT_S == m.FV_S + m.qF*m.FF_S
m.const_S19 = Constraint(rule=const_S19_rule)

def const_S20_rule(m):
    return m.F_comp == m.FMK + m.FT_S
m.const_S20 = Constraint(rule=const_S20_rule)

def const_S11_rule(m):
    # This equation specifies composition in "condenser," after 
    # mixing with make-up CO2 stream
    return m.FT_S*m.yS_CO2[4] + m.FMK == m.F_comp*m.yS_CO2[5]
m.const_S11 = Constraint(rule=const_S11_rule)

#def const_S12_rule(m): 
    # This equation specifies FSS var
#    return m.FVS+m.qF*m.FFS+m.FMK == m.FLS+m.FSS
#m.const_S12 = Constraint(rule=const_S12_rule)

def const_S12_rule(m):
    # Flow in the compressor equals combined reflux and distillate flow rates
    return m.F_comp == m.FL_S + m.FD_S
m.const_S12 = Constraint(rule=const_S12_rule)

# reboiler bottoms equation:
def const_S13_rule(m):
    return m.FB_S == m.FL_S + (1-m.qF)*m.FF_S - m.FV_S
m.const_S13 = Constraint(rule=const_S13_rule)

# VLE:
def const_S14_rule(m):
    return m.KS_CO2 == 1.3
m.const_S14 = Constraint(rule=const_S14_rule)

def const_S18_rule(m):
    return m.KS_IPA == 0.0213
m.const_S18 = Constraint(rule=const_S18_rule)

def const_S15_rule(m,i):
    return m.KS_IPA*m.EMVS*(1-m.xS_CO2[i]) == \
            (1-m.yS_CO2[i]) - (1-m.yS_CO2[i-1])*(1-m.EMVS)
m.const_S15 = Constraint(m.set_S_eq,rule=const_S15_rule)

# assume reboiler reaches full vapor liquid equilibrium:
def const_S16_rule(m):
    return m.KS_IPA*(1-m.xS_CO2[1]) == (1-m.yS_CO2[1])
m.const_S16 = Constraint(rule=const_S16_rule)

# compositions are the same across phases in the condenser
# (because they are different 'branches' of the same stream)
def const_S17_rule(m): 
    return m.yS_CO2[m.NTS.value] == m.xS_CO2[m.NTS.value]
m.const_S17 = Constraint(rule=const_S17_rule)

# Condenser:
def const_cond1_rule(m):
    return m.P1_cond == m.P2_F
m.const_cond1 = Constraint(rule=const_cond1_rule)
    
def const_cond2_rule(m):
    return m.T1_cond == m.TF
m.const_cond2 = Constraint(rule=const_cond2_rule)

def const_cond3_rule(m):
    return m.V1_cond == (1/m.P1_cond)*m.F_comp*0.08206*m.T1_cond
m.const_cond3 = Constraint(rule=const_cond3_rule)

def const_cond4_rule(m):
    return m.P2_cond == 100.0
m.const_cond4 = Constraint(rule=const_cond4_rule)

def const_cond5_rule(m):
    return m.P1_cond*m.V1_cond**1.4 == m.P2_cond*m.V2_cond**1.4
m.const_cond5 = Constraint(rule=const_cond5_rule)

def const_cond6_rule(m):
    return m.P2_cond*m.V2_cond == m.F_comp*0.08206*m.T2_cond
m.const_cond6 = Constraint(rule=const_cond6_rule)


# Reboiler
# reboiler already almost specified, just need energy balance
# need to:
#   enforce (mixture) heat capacity as a function of temp (actually not necessary for st.st.)
#   maybe enforce latent heat (of vaporization) of liquid from stripper (x2/L2?) (as fcn of Temp?)
#   calculate reboiler temperature
#   calculate temperature of vapor (solvent) leaving reboiler (with steadied en. bal. constraint)
#   ? calculate density of vapor CO2 solvent through reboiler ? (Ramchandran uses average)
#   ^ again not necessary for steady state

def const_B1_rule(m):
    return m.Pvap_B_IPA*101.325 == 10**(6.866 - 1360/(197.6+m.TB_sh-273.15))
m.const_B1 = Constraint(rule=const_B1_rule)

def const_B2_rule(m):
    return m.P_B == 40.0
m.const_B2 = Constraint(rule=const_B2_rule)

def const_B3_rule(m):
    return m.Pvap_B_IPA*(1-m.xS_CO2[1]) == m.P_B*(1-m.yS_CO2[1])
m.const_B3 = Constraint(rule=const_B3_rule)

def const_B4_rule(m):
    return m.dH_IPA == 44000 # kJ/kmol
m.const_B4 = Constraint(rule=const_B4_rule)

def const_B5_rule(m):
    return m.dH_CO2 == 15326 # kJ/kmol
m.const_B5 = Constraint(rule=const_B5_rule)

def const_B6_rule(m):
    return m.dH_B == m.yS_CO2[1]*m.dH_CO2 + (1-m.yS_CO2[1])*m.dH_IPA
m.const_B6 = Constraint(rule=const_B6_rule)

def const_B7_rule(m):
    return 0 == m.FV_S*m.dH_B - m.U_ov_B*m.A_B*(m.TB_tb_ave - m.TB_sh)
m.const_B7 = Constraint(rule=const_B7_rule)

def const_B8_rule(m):
    return m.CpB_CO2 == 40.0
m.const_B8 = Constraint(rule=const_B8_rule)

def const_B9_rule(m):
    return m.TB_tb_in == m.T2_cond
m.const_B9 = Constraint(rule=const_B9_rule)

def const_B10_rule(m):
    return (m.TB_tb_in - m.TB_tb_out)*m.CpB_CO2*m.F_B == m.FV_S*m.dH_B
m.const_B10 = Constraint(rule=const_B10_rule)

def const_B11_rule(m):
    return 2*m.TB_tb_ave - m.TB_tb_in == m.TB_tb_out
m.const_B11 = Constraint(rule=const_B11_rule)

def const_B12_rule(m):
    return m.F_B_bypass == m.FD_S - m.F_B
m.const_B12 = Constraint(rule=const_B12_rule)

# Cooler

def const_C1_rule(m):
    return 2*m.TC_tb_ave == m.TC_tb_out + m.TC_tb_in
m.const_C1 = Constraint(rule=const_C1_rule)

def const_C2_rule(m):
    return m.F_C == m.F_B + m.F_B_bypass
m.const_C2 = Constraint(rule=const_C2_rule)

def const_C3_rule(m):
    return m.TC_tb_in*m.F_C == m.TB_tb_out*m.F_B + m.T2_cond*m.F_B_bypass
m.const_C3 = Constraint(rule=const_C3_rule)

#def const_C4_rule(m):
#    return m.TC_cool == 306 
#m.const_C4 = Constraint(rule=const_C4_rule)

def const_C5_rule(m):
    return m.CpC_CO2 == m.CpB_CO2
m.const_C5 = Constraint(rule=const_C5_rule)

# assuming heat capacity of stream is that of CO2
def const_C6_rule(m):
    return 0 == m.F_C*m.CpC_CO2*(m.TC_tb_in-m.TC_tb_out) - m.U_ov_C*m.A_C*(m.TC_tb_ave-m.TC_sh_ave)
m.const_C6 = Constraint(rule=const_C6_rule)

def const_C7_rule(m):
    return m.CpC_H2O == 36.9
m.const_C7 = Constraint(rule=const_C7_rule)

def const_C8_rule(m):
    return m.TC_sh_in == 293
m.const_C8 = Constraint(rule=const_C8_rule)

def const_C9_rule(m):
    return 2*m.TC_sh_ave == m.TC_sh_in + m.TC_sh_out
m.const_C9 = Constraint(rule=const_C9_rule)

#def const_C10_rule(m):
#    return m.F_C_H2O == 22.25
#    #return m.F_C_H2O == 99.0 
#m.const_C10 = Constraint(rule=const_C10_rule)

def const_C11_rule(m):
    return 0 == m.F_C_H2O*m.CpC_H2O*(m.TC_sh_out-m.TC_sh_in) - m.U_ov_C*m.A_C*(m.TC_tb_ave-m.TC_sh_ave)
m.const_C11 = Constraint(rule=const_C11_rule)

def const_C12_rule(m):
    return m.yC_CO2 == m.xS_CO2[m.NTS]
m.const_C12 = Constraint(rule=const_C12_rule)

def const_C13_rule(m):
    return m.yC_IPA == 1 - m.yC_CO2
m.const_C13 = Constraint(rule=const_C13_rule)

# connectivity
def const_C14_rule(m):
    return m.yC_IPA == m.yE_S_IPA
m.const_C14 = Constraint(rule=const_C14_rule)

#def const_C15_rule(m):
#    return m.F_C == m.FSE
    #return m.FSE == m.FMK + m.FT_S - m.FD_S
    #return m.FMK + m.FV_S == m.FL_S
    #return m.FMK == m.FL_S-m.FV_S+(1-m.qF)*m.FF_S
#m.const_C15 = Constraint(rule=const_C15_rule)


###################
### objective #####
###################

#m.obj_ss = Objective(sense=minimize,expr=(m.TC_tb_out-313.0)**2)
m.obj_ss = Objective(sense=minimize,expr=(m.F_C-m.FSE)**2+(m.TC_tb_out-m.TE_sp)**2)
#m.obj_ss = Objective(sense=minimize,expr=0.0)



