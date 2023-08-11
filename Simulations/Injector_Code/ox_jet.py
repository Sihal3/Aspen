import numpy as np
from jet import Jet
import CoolProp.CoolProp as CP
from dimensions import ureg

class Ox_Jet(Jet):
    def __init__(self, id_num, package):
        super().__init__(package)
        self.package = package
        
        self.dia.append(package['ox_dia_u'].magnitude)
        self.id_num = id_num
        self.ox_temp = package['ox_temp_u'].magnitude
        self.oxidizer = package['oxidizer']
        self.inj_pressure = package['inj_pressure_u'].magnitude
        
    def add_density(self):
        den = CP.PropsSI('D', 'T', self.temp[-1], 'P', self.package['chamber_pressure_u'].magnitude, self.oxidizer)
        self.density.append(den)
    
    
    def add_phase(self):
        phase = CP.PhaseSI("P", self.inj_pressure, "T", self.ox_temp, self.oxidizer)
        
        if phase == "liquid":
            self.phase.append(1)
        else:
            self.phase.append(0) # assuming it can either be gas or liq

    
    def set_origins(self):
        # sets the initial positions for the jets 
        cup_ring_rad_u = self.cup_ring_dia_u / 2
        if self.id_num == 0:
            pos_x = 0
            pos_y = 0
        else:
            pos_x = cup_ring_rad_u * np.cos(self.id_num * (2 * np.pi) / (self.num_ox_core - 1)) 
            pos_y = cup_ring_rad_u * np.sin(self.id_num * (2 * np.pi) / (self.num_ox_core - 1)) 
        
        if pos_x != 0:
            pos_x = pos_x.magnitude

        if pos_y != 0:
            pos_y = pos_y.magnitude

        pos_z = - self.cup_inset_u
        self.pos_vectors.append(np.array([pos_x, pos_y, pos_z.magnitude]))
    
    def set_initial_velocity(self):
        # takes in a scalar velocity and makes it into a vector assuming it is parallel to chamber
        vel = self.package['ox_vel_u']
        self.vel_vectors.append(np.array([0, 0, vel.magnitude]))
    
    def calc_expansion(self):
        #TODO Update the diameter for the jet based on its position
        pass