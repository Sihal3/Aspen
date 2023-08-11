import matplotlib.pyplot as plt
import numpy as np
from dimensions import ureg
from fu_jet import Fu_Jet

# Film Jet is a type of Fu Jet that is used exclusively for 
#   film cooling the wall of the combustion chamber. 
# It extends the Fu_Jet class and adds in the property of 
#   phase for the film cooling jet. This is a necessary 
#   addition since film cooling utilizes liquid and gas
#   for creating a protective layer isolating the walls
#   from the combustion in the center of the chamber.
# Increased attention also needs to be paid to the heat
#   fluxes leaving the combustion area, since that will 
#   determine the effectiveness of the film cooling.

# Film cooling jets will always run parallel with the chamber axis

class Film_Jet(Fu_Jet):
    def __init__(self, id_num, package):
        super().__init__(id_num, package)
        self.film_ring_dia_u = package['film_ring_dia_u']
        self.num_film = package['num_film']
        self.film_vel_u = package['film_vel_u']
        
    def set_origins(self):
        pos_x = (self.film_ring_dia_u / 2) * np.cos(self.id_num * (2 * np.pi / self.num_film))
        pos_y = (self.film_ring_dia_u / 2) * np.sin(self.id_num * (2 * np.pi / self.num_film))
        pos_z = 0 * ureg('meter')
        
        self.pos_vectors.append(np.array([pos_x.magnitude, pos_y.magnitude, pos_z.magnitude]))
    
    def set_initial_velocity(self):
        self.vel_vectors.append(np.array([0, 0, self.film_vel_u.magnitude]))
