from jet import Jet
import numpy as np


class Fu_Jet(Jet):
    def __init__(self, id_num, package):
        super().__init__(package)
        self.package = package
        self.fu_angle_u = package['fu_angle_u']
        self.id_num = id_num

    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm

    def set_origins(self):
        
        if self.id_num < self.num_fu_per:

            pos_x = ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.cos(self.id_num * (2 * np.pi / self.num_fu_per))
            pos_y = ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.sin(self.id_num * (2 * np.pi / self.num_fu_per))
            
        
        else:
            ox_location = np.floor(self.id_num / self.num_fu_per)
            fu_location = self.id_num % self.num_fu_per
            
            pos_x = (self.cup_ring_dia_u / 2) * np.cos(ox_location * (2 * np.pi / (self.num_ox_core - 1))) + ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.cos(fu_location * (2 * np.pi / self.num_fu_per))
            pos_y = (self.cup_ring_dia_u / 2) * np.sin(ox_location * (2 * np.pi / (self.num_ox_core - 1))) + ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.sin(fu_location * (2 * np.pi / self.num_fu_per))
                
        pos_z = -(self.cup_inset_u / 2)        
                
        self.pos_vectors.append(np.array([pos_x.magnitude, pos_y.magnitude, pos_z.magnitude]))


    
    def set_initial_velocity(self):
        # set vector for velocity
        fu_location = self.id_num % self.num_fu_per
        vel = self.package['fu_vel_u'].magnitude
        ang = (self.fu_angle_u - np.radians(270)).magnitude
        psi = (((2 * np.pi) / self.num_fu_per) * fu_location)
        direction_vector = np.array([np.cos(psi) * np.cos(ang), np.sin(psi) * np.cos(ang), np.sin(ang)])
        unit_vel = vel * self.normalize(direction_vector)
        self.vel_vectors.append(unit_vel)
        
