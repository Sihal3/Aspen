import numpy as np
from jet import Jet

from dimensions import ureg

class Ox_Jet(Jet):
    def __init__(self, id_num, package):
        super().__init__(package)
        self.package = package
        
        self.dia.append(package['ox_dia_u'])
        self.id_num = id_num
        
    
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
    