from pint import UnitRegistry
import numpy as np
from CoolProp.CoolProp import PropsSI

ureg = UnitRegistry()


def ox_massflow(B, mDot):
    return (B * mDot) / (B + 1)

def fu_massflow(B, mDot):
    return mDot / (B + 1)

def incompressible_massflow(Cd, p, P_1, P_2, d, D):
    C = (1 /  (Cd - (d / D) ** 4)) ** 0.5
    A2 = (d / 2) ** 2 * np.pi
    massflow = C * A2 * (2 * p * (P_1 - P_2)) ** 0.5
    return massflow

def incompressible_orifice(Cd, p, mDot, P_1, P_2):
    area = mDot / (Cd * ((2 * p * (P_1 - P_2)) ** 0.5))
    diameter = 2 * (area / np.pi) ** 0.5
    return diameter


def fluid_vel(rho, dia, mDot):
    area = (dia / 2) ** 2 * np.pi
    return mDot / (rho * area)

def within_tol(dia, drill_dia, tolerance):
    # checks to see if the orifice diameter is close enough to the diameter of the drill
    upper_dia = (dia + tolerance).magnitude
    lower_dia = (dia - tolerance).magnitude
    if upper_dia > drill_dia.magnitude and lower_dia < drill_dia.magnitude:
        return True
    return False


def calc(inputs):
    # Unpack the input dictionary, allows expansion without difficulty
    chamber_pressure_u = inputs['chamber_pressure_u']
    B = inputs['o_f_ratio']
    mDot_u = inputs['mDot_u']
    ox_temp_u = inputs['ox_temp_u']
    fu_temp_u = inputs['fu_temp_u']
    fuel = inputs['fuel']
    oxidizer = inputs['oxidizer']
    film_cooling_percent = inputs['film_cooling_percent']
    cup_inset_u = inputs['cup_inset_u']
    cup_ring_dia_u = inputs['cup_ring_dia_u']
    chamber_dia_u = inputs['chamber_dia_u']
    chamber_len_u = inputs['chamber_len_u']
    film_wall_distance_u = inputs['film_wall_distance_u']
    num_ox_core = inputs['num_ox_core']
    fu_angle_u = inputs['fu_angle_u']
    res = inputs['res']
    
    # convert the film_wall_dist to a diameter for the film cooling ring
    film_cool_ring_u = chamber_dia_u - (2 * film_wall_distance_u)
    
    # list of drill diameters to reference
    drill_dia_u = {'1/64' : 0.396875 * ureg('mm'), '1/32' : 0.79375 * ureg('mm'), '1/16' : 1.5875 * ureg('mm')}
    
    pressure_drop_injector = .2 # comes from the standard 20% pressure drop through an injector
    inj_pressure_u = chamber_pressure_u * (pressure_drop_injector + 1)
    
    fu_mDot_u = fu_massflow(B, mDot_u)
    ox_mDot_u = ox_massflow(B, mDot_u)
    film_mDot_u = fu_mDot_u * film_cooling_percent
    
    Cd = .6 # set the Cd to something reasonable for the orifices  
    
    ox_density_u = PropsSI('D', 'T', ox_temp_u.magnitude, 'P', inj_pressure_u.magnitude, oxidizer) * ureg('kg/m^3')
    fu_density_u = PropsSI('D', 'T', fu_temp_u.magnitude, 'P', inj_pressure_u.magnitude, fuel) * ureg('kg/m^3')
    
    # section that determines the best number of fu orifices such that the diameter
    # is close to a drill bit diameter
    finding_dia = True
    fu_element_total = 0
    fu_element_dia_u = -1
    
    while(finding_dia):
        fu_element_total += 1
        fu_element_mDot_u = fu_mDot_u / fu_element_total # mDot through the current number of elements
        
        target_drill_u = drill_dia_u['1/32'].to(ureg.m)
        tolerance_u = 0.1 * ureg('mm')
        tolerance_u.ito(ureg.m)
        fu_element_dia_u = incompressible_orifice(Cd, fu_density_u, fu_element_mDot_u, inj_pressure_u, chamber_pressure_u)
        finding_dia = not within_tol(fu_element_dia_u, target_drill_u, tolerance_u)

    finding_dia = True
    film_dia_u = -1
    num_film = 0
    while(finding_dia):
        num_film += 1
        film_element_mDot_u = film_mDot_u / num_film
        
        film_target_drill_u = drill_dia_u['1/32'].to(ureg.m)
        tolerance_u = 0.1 * ureg('mm')
        tolerance_u.ito(ureg.m)
        film_dia_u = incompressible_orifice(Cd, fu_density_u, film_element_mDot_u, inj_pressure_u, chamber_pressure_u)
        finding_dia = not within_tol(film_dia_u, film_target_drill_u, tolerance_u)


    fu_elements_per = np.floor(fu_element_total / (num_ox_core))
    fu_elements_center = (fu_element_total - fu_elements_per * num_ox_core) + fu_elements_per

    ox_element_mDot_u = ox_mDot_u / num_ox_core
    ox_element_dia_u = incompressible_orifice(Cd, ox_density_u, ox_element_mDot_u, inj_pressure_u, chamber_pressure_u)


    # the cup is just going to be 1.x times bigger than the ox orifice for now
    cup_dia_u = ox_element_dia_u * 1.25

    # jet velocities
    ox_element_vel_u = fluid_vel(ox_density_u, ox_element_dia_u, ox_element_mDot_u)
    fu_element_vel_u = fluid_vel(fu_density_u, fu_element_dia_u, fu_element_mDot_u)
    film_vel_u = fluid_vel(fu_density_u, film_dia_u, film_element_mDot_u)
    
    
    
    
    
    # packaging    
    ans = { 'num_fu_per' : fu_elements_per, 
            'num_fu_center' : fu_elements_center, 
            'fu_vel_u' : fu_element_vel_u.to_base_units(), 
            'ox_vel_u' : ox_element_vel_u.to_base_units(),
            'num_ox_core' : num_ox_core,
            'fu_dia_u' : fu_element_dia_u.to_base_units(),
            'ox_dia_u' : ox_element_dia_u.to_base_units(),
            'fu_density_u' : fu_density_u.to_base_units(),
            'ox_density_u' : ox_density_u.to_base_units(),
            'fuel' : fuel,
            'oxidizer' : oxidizer,
            'chamber_pressure_u' : chamber_pressure_u.to_base_units(),
            'fu_temp_u' : fu_temp_u.to_base_units(),
            'ox_temp_u' : ox_temp_u.to_base_units(),
            'fu_element_total' : fu_element_total,
            'fu_mDot_u' : fu_mDot_u,
            'ox_mDot_u' : ox_mDot_u,
            'o_f_ratio' : B,
            'film_dia_u' : film_dia_u.to_base_units(),
            'num_film' : num_film,
            'film_mDot_u' : film_mDot_u,
            'film_vel_u' : film_vel_u.to_base_units(),
            'cup_inset_u' : cup_inset_u,
            'cup_ring_dia_u' : cup_ring_dia_u,
            'chamber_dia_u' : chamber_dia_u,
            'chamber_len_u' : chamber_len_u,
            'film_ring_dia_u' : film_cool_ring_u,
            'fu_angle_u' : fu_angle_u,
            'cup_dia_u' : cup_dia_u,
            'res' : res}
    
    return ans

def disp(calc_data):
    # unpacking
    fuel = calc_data['fuel']
    oxidizer = calc_data['oxidizer']
    fu_vel_u = calc_data['fu_vel_u']
    ox_vel_u = calc_data['ox_vel_u']
    num_ox_core = calc_data['num_ox_core']
    num_fu_per = calc_data['num_fu_per']
    num_fu_center = calc_data['num_fu_center']
    fu_dia_u = calc_data['fu_dia_u']
    ox_dia_u = calc_data['ox_dia_u']
    fu_density_u = calc_data['fu_density_u']
    ox_density_u = calc_data['ox_density_u']
    chamber_pressure_u = calc_data['chamber_pressure_u']
    fu_temp_u = calc_data['fu_temp_u']
    ox_temp_u = calc_data['ox_temp_u']
    fu_elements_total = calc_data['fu_element_total']
    fu_mDot_u = calc_data['fu_mDot_u']
    ox_mDot_u = calc_data['ox_mDot_u']
    o_f_ratio = calc_data['o_f_ratio']
    film_dia_u = calc_data['film_dia_u']
    num_film = calc_data['num_film']
    film_mDot_u = calc_data['film_mDot_u']
    film_vel_u = calc_data['film_vel_u']
    cup_ring_dia_u = calc_data['cup_ring_dia_u']
    cup_inset_u = calc_data['cup_inset_u']
    
    
    print("\n\n\n")
    print(f"With {fuel} & {oxidizer} propellants, and under the following conditions...\n")
    print(f"    Chamber pressure of {chamber_pressure_u.to(ureg.psi):.4}")
    print(f"    O/F ratio of {o_f_ratio}")
    print(f"        {fuel} temperature : {fu_temp_u}")
    print(f"        {fuel} density : {fu_density_u:.4}")
    print(f"        {oxidizer} temperature : {ox_temp_u}")
    print(f"        {oxidizer} density : {ox_density_u:.4}")
    print()
    print(f"Results in the following injector characteristics...")
    print()
    print(f"    Total {fuel} massflow : {fu_mDot_u:.4}")
    print(f"    Total {oxidizer} massflow : {ox_mDot_u:.4}")
    print(f"    Total film cooling massflow : {film_mDot_u:.4}")
    print()
    print(f"    {fuel} orifice diameter : {fu_dia_u.to(ureg.inch):.3}")
    print(f"    Total {fuel} orifices : {fu_elements_total}")
    print(f"    Number of {fuel} orifices on outer ring : {num_fu_per}")
    print(f"    Number of {fuel} orifices in center core : {num_fu_center}")
    print(f"    {fuel} jet velocity : {fu_vel_u:.4}")
    print()
    print(f"    Diameter of film cooling orifices : {film_dia_u.to(ureg('inch')):.3}")
    print(f"    Total film cooling orifices : {num_film}")
    print(f"    Film cooling jet velocity : {film_vel_u:.4}")
    print()
    print(f"    {oxidizer} orifice diameter : {ox_dia_u.to(ureg.inch):.3}")
    print(f"    Number of {oxidizer} orifices : {num_ox_core}")
    print(f"    {oxidizer} jet velocity : {ox_vel_u:.4}")
    print("\n\n\n")
