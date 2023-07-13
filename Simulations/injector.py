import math
from pint import UnitRegistry
from CoolProp.CoolProp import PropsSI


ureg = UnitRegistry()


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
    return mDot / (Cd * (2 * g * p * dP)**.5)
    

def incompressible_massflow(Cd, p, P_1, P_2, d, D):
    C = (1 /  (Cd - (d / D) ** 4)) ** 0.5
    A2 = (d / 2) ** 2 * math.pi
    massflow = C * A2 * (2 * p * (P_1 - P_2)) ** 0.5
    return massflow

def incompressible_orifice(Cd, p, mDot, P_1, P_2):
    area = mDot / (Cd * ((2 * p * (P_1 - P_2)) ** 0.5 ))

    diameter = 2 * (area / math.pi) ** 0.5
    return diameter





def ox_injector_rad(Cd, p, dP, B, mDot, g):
    den = Cd * math.pi * (2 * g * p * dP)**.5 * (B + 1)
    num = mDot * B
    return (num / den)**.5

def fu_injector_rad(Cd, p, dP, B, mDot, g, Cd_o, p_o, dP_o, w_oi):
    term1 = mDot / (Cd * (2 * g * p * dP)**.5 * (B + 1))
    term2 = (ox_injector_rad(Cd_o, p_o, dP_o, B, mDot, g) + w_oi) ** 2 * math.pi
    return ((term1 + term2) / math.pi)**.5 

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

massflow = 0.8 * ureg('kilogram / second')
B = 3 
chamber_pressure = 180 * ureg('psi')
chamber_pressure.ito('pascal')
Cd = .7
film_cooling_percentage = .10





nitrous_mass_flow = ox_massflow(B, massflow)
ethanol_mass_flow = fu_massflow(B, massflow) * (1 + film_cooling_percentage)


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





# impinging SECTION

def inj_func(mDot, B, chamber_pressure, ign_chamber_press):
    print()
    drill_dia = {'1/64' : 0.396875 * ureg('mm'), '1/32' : 0.79375 * ureg('mm'), '1/16' : 1.5875 * ureg('mm')}
    
    ox_cores = 6
    fu_inj_per = 5
    film_cooling_orifices = 12
    
    fu_temp = 290 * ureg('kelvin')
    ox_temp = 244.15 * ureg('kelvin')

    fu_inj_pressure = chamber_pressure * 1.2
    
    fu_mDot = fu_massflow(B, mDot)
    print(f"Ethanol massflow rate: {fu_mDot:.4f}")
    
    ox_mDot = ox_massflow(B, mDot)
    print(f"Nitrous massflow rate: {ox_mDot:.4f}")
    
    film_cooling_mDot = fu_mDot * film_cooling_percentage
    print(f'Film cooling massflow: {film_cooling_mDot:.3f}')
    
    
    ox_density = PropsSI('D', 'T', ox_temp.magnitude, 'P', ign_chamber_press.magnitude, 'N2O') * ureg('kg/m^3')
    fu_density = PropsSI('D', 'T', fu_temp.magnitude, 'P', fu_inj_pressure.magnitude, "Ethanol") * ureg('kg/m^3')
    
    fu_diameter = incompressible_orifice(Cd, fu_density, fu_mDot, fu_inj_pressure, chamber_pressure)
    ox_diameter = incompressible_orifice(Cd, ox_density, ox_mDot, fu_inj_pressure, chamber_pressure)
    print()
    print(f"Ethanol density: {fu_density:.3f}")
    print(f"Nitrous density: {ox_density:.3f}")
    print()

    element_fu_mDot = fu_mDot / (ox_cores * fu_inj_per)
    element_ox_mDot = ox_mDot / ox_cores
    
    print(f'Ethanol massflow thru element: {element_fu_mDot:.4f}')
    print(f"Nitrous Massflow thru element: {element_ox_mDot:.4f}")
    print()
    element_fu_dia = incompressible_orifice(Cd, fu_density, element_fu_mDot, fu_inj_pressure, chamber_pressure).to_base_units()
    element_ox_dia = incompressible_orifice(Cd, ox_density, element_ox_mDot, fu_inj_pressure, chamber_pressure).to_base_units()
    print(f"Ethanol diameter with {fu_inj_per} elements per ox core: {element_fu_dia.to(ureg.inch):.3f}")
    print(f"Nitrous diameter with {ox_cores} elements: {element_ox_dia.to(ureg.inch):.3f}")
    
    film_cooling_dia = incompressible_orifice(Cd, fu_density, film_cooling_mDot / film_cooling_orifices, fu_inj_pressure, chamber_pressure)
    print(f"Film cooling diameter for {film_cooling_orifices} elements: {film_cooling_dia.to_base_units().to(ureg.inch):.3f}")
    
    print()
    # for i in range(1,29):
    #     element_fu_mDot = fu_mDot / i
    #     element_ox_mDot = ox_mDot / i
    #     print(f"Ethanol Massflow thru element: {element_fu_mDot:.4f}")
    #     print(f"Nitrous Massflow thru element: {element_ox_mDot:.4f}")
    #     element_fu_dia = incompressible_orifice(Cd, fu_density, element_fu_mDot, fu_inj_pressure, chamber_pressure).to_base_units()
    #     element_ox_dia = incompressible_orifice(Cd, ox_density, element_ox_mDot, fu_inj_pressure, chamber_pressure).to_base_units()
    #     print(f"Ethanol diameter with {i} elements: {element_fu_dia.to(ureg.mm):.3f}")
    #     print(f"Nitrous diameter with {i} elements: {element_ox_dia.to(ureg.mm):.3f}")    
    #     print()
    

inj_func(massflow, B, chamber_pressure,  chamber_pressure)





