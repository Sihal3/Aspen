from jet import Jet
import numpy as np


class Fu_Jet(Jet):
    def __init__(self, id_num, num_ox_core, num_fu_per, cup_inset, cup_ring_dia, cup_dia, fu_dia, fu_angle):
        super().__init__(num_ox_core, num_fu_per, cup_inset, cup_ring_dia, cup_dia)
        self.fu_dia = fu_dia
        self.fu_angle = fu_angle
        self.id_num = id_num
        self.dia.append(self.fu_dia)



    def set_origins(self):
        
        if self.id_num < self.num_fu_per:

            pos_x = ((self.cup_inset / (np.tan(self.fu_angle) * 2)) + (self.cup_dia / 2)) * np.cos(self.id_num * (2 * np.pi / self.num_fu_per))
            pos_y = ((self.cup_inset / (np.tan(self.fu_angle) * 2)) + (self.cup_dia / 2)) * np.sin(self.id_num * (2 * np.pi / self.num_fu_per))
            
        
        else:
            ox_location = np.floor(self.id_num / self.num_fu_per)
            fu_location = self.id_num % self.num_fu_per
            
            pos_x = (self.cup_ring_dia / 2) * np.cos(ox_location * (2 * np.pi / (self.num_ox_core - 1))) + ((self.cup_inset / (np.tan(self.fu_angle) * 2)) + (self.cup_dia / 2)) * np.cos(fu_location * (2 * np.pi / self.num_fu_per))
            pos_y = (self.cup_ring_dia / 2) * np.sin(ox_location * (2 * np.pi / (self.num_ox_core - 1))) + ((self.cup_inset / (np.tan(self.fu_angle) * 2)) + (self.cup_dia / 2)) * np.sin(fu_location * (2 * np.pi / self.num_fu_per))
                
        pos_z = -(self.cup_inset / 2)        
                
        self.pos_vectors.append(np.array([pos_x, pos_y, pos_z]))


    
    def set_initial_velocity(self, vel):
        # set vector for velocity
        fu_location = self.id_num % self.num_fu_per
        
        ang = self.fu_angle - np.radians(270) 
        psi = ((2 * np.pi) / self.num_fu_per) * fu_location
        direction_vector = np.array([np.cos(psi) * np.cos(ang), np.sin(psi) * np.cos(ang), np.sin(ang)])
        unit_vel = vel * normalize(direction_vector)
        self.vel_vectors.append(unit_vel)
        
