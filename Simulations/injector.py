import math
import pint


# -- Definitions -- 
# mDot = mass flow
# p = density (rho)
# dp = delta Pressure
# g = gravitational constant
# B = o/f ratio
# Cd = coefficient
# w_.. = width
# A = cross-sectional area


def orifice_equation(Cd, p, dP, mDot, g):
    return mDot / (Cd * math.sqrt(2 * g * p * dP))
    
    
def ox_injector_rad(Cd, p, dP, B, mDot, g):
    den = Cd * math.pi * math.sqrt(2 * g * p * dP) * (B + 1)
    num = mDot * B
    return math.sqrt(num / den)

def fu_injector_rad(mDot, Cd, g, p, dP, B, Cd_o, p_o, dP_o, w_oi):
    term1 = mDot / (Cd * math.sqrt(2 * g * p * dP) * (B + 1))
    term2 = (ox_injector_rad(Cd_o, p_o, dP_o, B, mDot, g) + w_oi) ** 2 * math.pi
    return math.sqrt((term1 + term2) / math.pi) 

def recess_len(w_oi, Cd, p, dP, B, mDot, g):
    return 4 * (w_oi + ox_injector_rad(Cd, p, dP, B, mDot, g))


def fu_massflow(B, mDot):
    return mDot / (B + 1)

def ox_massflow(B, mDot):
    return (B * mDot) / (B + 1)


def line_velocity(mDot, p, A):
    return mDot / (p * A)






# OVERARCHING SECTION
# section to input all of the necessary inputs that cascade through the 
# whole system

massflow = 0.69 
B = 3
chamber_pressure = 180
Cd = .9
film_cooling_percentage = .10
nitrous_density = 200
ethanol_density = 200





nitrous_mass_flow = ox_massflow(B, massflow)
ethanol_mass_flow = fu_massflow(B, massflow) * (1 + film_cooling_percentage)


# IGNECTOR SECTION
# section outlining the details for the inner core of the coaxial injector
# which will supply the nitrous to the engine, while also serving as the 
# ignition source for startup.




# outputs
ox_pressure = -1
fu_pressure = -1
ox_fitting_num = -1 # output the fitting number for the fitting leading nitrous into the inner core

fu_ign_massflow = .003

ignector_pressure = 1.2 * chamber_pressure # 20% more than chamber pressure, as per the usual 
fu_pressure = 1.2 * ignector_pressure 

# determine the best combination for orifice size and pressure
fitting_dia = {4 : .25, 6 : .375, 8 : .5, 10 : .675}


ox_dia = orifice_equation(Cd, nitrous_density, )
fu_dia = ((orifice_equation(Cd, ethanol_density, fu_pressure, fu_ign_massflow, 9.81) / 2) ** 2) / math.pi






# COAXIAL SECTION

def coax_func(mDot, B, chamber_pressure, ign_chamber_press):
    fu_inj_pressure = chamber_pressure * 1.2
    
    fu_injector_rad(mDot, Cd, 9.81, ethanol_density, fu_inj_pressure - chamber_pressure, B, Cd, nitrous_density, )