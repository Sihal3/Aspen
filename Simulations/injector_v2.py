import math
from pint import UnitRegistry
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt


# import the necessary files related to the simulation
import dimensions, jet, fu_jet, ox_jet, sim

# start pint to keep track of units throughout the simulation
from dimensions import ureg
Q_ = ureg.Quantity


# Simulation inputs
# enter values in imperial (lbs, psi, deg, etc)
# temp MUST be in kelvin because i am not that insane
massflow = 2.2
o_f_ratio = 1.2
chamber_pressure = 180
film_cooling_percentage = .15
fuel = 'Ethanol'
oxidizer = 'N2O'
fu_inj_angle = 50
fu_temp = 270
ox_temp = 260
cup_inset = .2
cup_ring_dia = 1.5 
chamber_dia = 3
chamber_len = 10 # length of the chamber
film_wall_distance = 0.1 # distance from the wall to a film cooling orifice
num_ox_core = 7     # set the number of oxidizer injectors because I am too lazy to make it automated    



# Adding units to the inputs and converting to SI
massflow_u = massflow * ureg('pound / second')
massflow_u.ito('kilogram / second')

chamber_pressure_u = chamber_pressure * ureg('psi')
chamber_pressure_u.ito('pascal')

fu_inj_angle_u = fu_inj_angle * ureg('degree')
fu_inj_angle_u.ito('radians')

fu_temp_u = fu_temp * ureg('kelvin')

ox_temp_u = ox_temp * ureg('kelvin')

cup_inset_u = cup_inset * ureg('inch')
cup_inset_u.ito('meter')

cup_ring_dia_u = cup_ring_dia * ureg('inch')
cup_ring_dia_u.ito('meter')

chamber_dia_u = chamber_dia * ureg('inch')
chamber_dia_u.ito('meter')

chamber_len_u = chamber_len * ureg('inch')
chamber_len_u.ito('meter')

film_wall_distance_u = film_wall_distance * ureg('inch')
film_wall_distance_u.ito('meter')

# Determining the massflow for each propellant
fu_massflow_u = dimensions.fu_massflow(o_f_ratio, massflow_u)
ox_massflow_u = dimensions.ox_massflow(o_f_ratio, massflow_u)
film_massflow_u = fu_massflow_u * (1 + film_cooling_percentage)



package = { 'mDot_u' : massflow_u, 
            'chamber_pressure_u' : chamber_pressure_u, 
            'o_f_ratio' : o_f_ratio,
            'fu_temp_u' : fu_temp_u,
            'ox_temp_u' : ox_temp_u,
            'fuel' : fuel,
            'oxidizer' : oxidizer,
            'film_cooling_percent' : film_cooling_percentage,
            'fu_angle_u' : fu_inj_angle_u,
            'cup_ring_dia_u' : cup_ring_dia_u,
            'cup_inset_u' : cup_inset_u,
            'chamber_dia_u' : chamber_dia_u,
            'chamber_len_u' : chamber_len_u,
            'film_wall_distance_u' : film_wall_distance_u,
            'num_ox_core' : num_ox_core}


out = dimensions.calc(package)
dimensions.disp(out)

simulation = sim.Sim()

simulation.create_scenario(out)
simulation.make_injector()
simulation.go()
simulation.disp()