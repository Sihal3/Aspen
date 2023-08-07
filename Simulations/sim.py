import numpy as np
from pint import UnitRegistry
import matplotlib.pyplot as plt

import dimensions
from dimensions import ureg

from ox_jet import Ox_Jet
from fu_jet import Fu_Jet



class Sim:
    def __init__(self):
        self.fu_list = []
        self.ox_list = []

        plt.rcParams['figure.figsize'] = (10,10)
        self.ax = plt.axes(projection='3d')
        self.ax.set_xlim(-10, 10) 
        self.ax.set_ylim(-10, 10)
        self.ax.set_zlim(-10, 10)
    

    def create_scenario(self, calc_data):
        self.fuel = calc_data['fuel']
        self.oxidizer = calc_data['oxidizer']
        self.fu_vel_u = calc_data['fu_vel_u']
        self.ox_vel_u = calc_data['ox_vel_u']
        self.num_ox_core = int(calc_data['num_ox_core'])
        self.num_fu_per = int(calc_data['num_fu_per'])
        self.num_fu_center = int(calc_data['num_fu_center'])
        self.fu_dia_u = calc_data['fu_dia_u']
        self.ox_dia_u = calc_data['ox_dia_u']
        self.fu_density_u = calc_data['fu_density_u']
        self.ox_density_u = calc_data['ox_density_u']
        self.chamber_pressure_u = calc_data['chamber_pressure_u']
        self.fu_temp_u = calc_data['fu_temp_u']
        self.ox_temp_u = calc_data['ox_temp_u']
        self.fu_elements_total = calc_data['fu_element_total']
        self.fu_mDot_u = calc_data['fu_mDot_u']
        self.ox_mDot_u = calc_data['ox_mDot_u']
        self.o_f_ratio = calc_data['o_f_ratio']
        self.film_dia_u = calc_data['film_dia_u']
        self.num_film = calc_data['num_film']
        self.film_mDot_u = calc_data['film_mDot_u']
        self.film_vel_u = calc_data['film_vel_u']
        self.chamber_dia_u = calc_data['chamber_dia_u']
        self.chamber_len_u = calc_data['chamber_len_u']
        self.film_ring_dia_u = calc_data['film_ring_dia_u']
        self.cup_ring_dia_u = calc_data['cup_ring_dia_u']
        self.cup_inset_u = calc_data['cup_inset_u']
        self.fu_angle_u = calc_data['fu_angle_u']
        self.cup_dia_u = calc_data['cup_dia_u']
        self.res = calc_data['res']
        
        
        for ox_id in range(self.num_ox_core):
            ox_jet = Ox_Jet(ox_id, calc_data)
            ox_jet.set_initial_velocity()
            ox_jet.set_origins()
            ox_jet.set_acc()
            
            self.ox_list.append(ox_jet)
        
        for fu_id in range(self.fu_elements_total):
            fu_jet = Fu_Jet(fu_id, calc_data)
            fu_jet.set_origins()
            fu_jet.set_initial_velocity()
            fu_jet.set_acc()
            
            self.fu_list.append(fu_jet)
            
    #TODO make the create scene function take in and set all the variables in the package
    #TODO and fill the list of the fu and ox jets 


    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm
    
    def plot_rotating_ellipse_3d(self, a, b, direction_vector, center_point, theta):
        # Normalize the direction vector to ensure it is a unit vector
        direction_vector = direction_vector / np.linalg.norm(direction_vector)
        
        # Convert theta from degrees to radians
        theta = np.radians(theta)
        
        # Generate points on the ellipse in 3D using parametric equations
        t = np.linspace(0, 2 * np.pi, 100)
        
        # Generate the normal vector to the plane of the ellipse
        n = np.cross(direction_vector, np.array([0, 0, 1]))
        n = n / np.linalg.norm(n)
        
        # Rotate the direction vector around its axis (theta)
        rotation_matrix = np.array([[np.cos(theta) + direction_vector[0]**2 * (1 - np.cos(theta)),
                                    direction_vector[0] * direction_vector[1] * (1 - np.cos(theta)) - direction_vector[2] * np.sin(theta),
                                    direction_vector[0] * direction_vector[2] * (1 - np.cos(theta)) + direction_vector[1] * np.sin(theta)],
                                    [direction_vector[0] * direction_vector[1] * (1 - np.cos(theta)) + direction_vector[2] * np.sin(theta),
                                    np.cos(theta) + direction_vector[1]**2 * (1 - np.cos(theta)),
                                    direction_vector[1] * direction_vector[2] * (1 - np.cos(theta)) - direction_vector[0] * np.sin(theta)],
                                    [direction_vector[0] * direction_vector[2] * (1 - np.cos(theta)) - direction_vector[1] * np.sin(theta),
                                    direction_vector[1] * direction_vector[2] * (1 - np.cos(theta)) + direction_vector[0] * np.sin(theta),
                                    np.cos(theta) + direction_vector[2]**2 * (1 - np.cos(theta))]])
        
        direction_vector = np.dot(rotation_matrix, direction_vector)
        
        # Generate the ellipse points in 3D space using the rotated direction vector and normal vector
        x = center_point[0] + a * (np.cos(t) * direction_vector[0] + np.sin(t) * n[0])
        y = center_point[1] + b * (np.cos(t) * direction_vector[1] + np.sin(t) * n[1])
        z = center_point[2] + a * direction_vector[2] * np.cos(t) + b * n[2] * np.sin(t)
        
        # Plot the ellipse
        self.ax.plot(x, y, z, color='b', alpha=0.6)

    def plot_rotating_circle_3d(self, radius, direction_vector, center, color= 'black'):
        # radius is int
        # direction_vector is an np array
        # center is an np array representing a point
        
        # Normalize the direction vector to ensure it is a unit vector
        direction_vector = direction_vector / np.linalg.norm(direction_vector)
        
        # Generate points on the circle in 3D using parametric equations
        t = np.linspace(0, 2 * np.pi, 100)
        # Generate the normal vector to the plane of the circle
        n = np.cross(direction_vector, np.array([0, 0, 1]))
        n = self.normalize(n)
        
        # Generate the circle points in 3D space using the direction vector and normal vector
        x = center[0] + radius * (np.cos(t) * direction_vector[0] + np.sin(t) * n[0])
        y = center[1] + radius * (np.cos(t) * direction_vector[1] + np.sin(t) * n[1])
        z = center[2] + radius * (np.cos(t) * direction_vector[2] + np.sin(t) * n[2])
        
        
        
        # Plot the circle
        plt.plot(x, y, z, color=color, alpha=0.6)
    
    def make_injector(self):
        # convert all of the length dimensions to inches since that is the units for the plot
        chamber_dia = self.chamber_dia_u.to(ureg('inch')).magnitude
        chamber_len = self.chamber_len_u.to(ureg('inch')).magnitude
        cup_ring_dia = self.cup_ring_dia_u.to(ureg('inch')).magnitude
        cup_inset = self.cup_inset_u.to(ureg('inch')).magnitude
        fu_dia = self.fu_dia_u.to(ureg('inch')).magnitude
        ox_dia = self.ox_dia_u.to(ureg('inch')).magnitude
        fu_angle = self.fu_angle_u.magnitude
        cup_dia = self.cup_dia_u.to(ureg('inch')).magnitude
        num_ox_core = self.num_ox_core
        num_fu_center = self.num_fu_center
        num_fu_per = self.num_fu_per
        film_ring_dia = self.film_ring_dia_u.to(ureg('inch')).magnitude
        num_film = self.num_film
        film_dia = self.film_dia_u.to(ureg('inch')).magnitude


        r = chamber_dia / 2
        
        theta = np.linspace(0, 2 * np.pi, 100)
        injector_face_xs = r * np.cos(theta)
        injector_face_ys = r * np.sin(theta)
        
        plt.plot(injector_face_xs, injector_face_ys, 0, color = 'black')
        plt.plot(injector_face_xs, injector_face_ys, chamber_len, color = 'black')
        cup_z = -cup_inset
        
        # center injector
        
            # outline flush with the injector
        center_cup_flush_dia = cup_dia + 2 * (cup_inset / np.tan(fu_angle))
        center_cup_flush_xs = (center_cup_flush_dia / 2) * np.cos(theta)
        center_cup_flush_ys = (center_cup_flush_dia / 2) * np.sin(theta)
        center_cup_flush_z = 0
        
        plt.plot(center_cup_flush_xs, center_cup_flush_ys, center_cup_flush_z, color = 'black')
        
            # cup inset
        center_cup_xs = (cup_dia / 2) * np.cos(theta)
        center_cup_ys = (cup_dia / 2) * np.sin(theta)
        center_cup_z = - cup_inset
        
        plt.plot(center_cup_xs, center_cup_ys, center_cup_z, color = 'black')
        
            # ox orifice
        center_ox_orf_xs = (ox_dia / 2) * np.cos(theta)
        center_ox_orf_ys = (ox_dia / 2) * np.sin(theta)
        center_ox_orf_z = - cup_inset
        
        plt.plot(center_ox_orf_xs, center_ox_orf_ys, center_ox_orf_z, color = 'black')
        
        
        for i in range(0, num_ox_core - 1):
            # plot the 3 circles representing the cups 
            
            # the offset for the cups that creates the cup pattern
            cup_ring_rad = cup_ring_dia / 2
            offset_x = cup_ring_rad * np.cos(i * (2 * np.pi) / (num_ox_core - 1))
            offset_y = cup_ring_rad * np.sin(i * (2 * np.pi) / (num_ox_core - 1))
            
            # the first one is the outline flush with the injector face
            
            cup_flush_dia = cup_dia + 2 * (cup_inset / np.tan(fu_angle))
            
            
            cup_flush_xs = (cup_flush_dia / 2) * np.cos(theta)
            cup_flush_ys = (cup_flush_dia / 2) * np.sin(theta)
            cup_flush_z = 0
            
            cup_flush_xs += offset_x 
            cup_flush_ys += offset_y
            
            plt.plot(cup_flush_xs, cup_flush_ys, cup_flush_z, color = 'black')
            
            # the second is the inset 
            cup_xs = (cup_dia / 2) * np.cos(theta)
            cup_ys = (cup_dia / 2) * np.sin(theta)
            cup_z = - cup_inset
            
            cup_xs += offset_x
            cup_ys += offset_y
            
            plt.plot(cup_xs, cup_ys, cup_z, color = 'black')
            
            # the third is the oxidizer orifice
            ox_orf_xs = (ox_dia / 2) * np.cos(theta)
            ox_orf_ys = (ox_dia / 2) * np.sin(theta)
            ox_orf_z = - cup_inset
            
            ox_orf_xs += offset_x
            ox_orf_ys += offset_y
            
            plt.plot(ox_orf_xs, ox_orf_ys, ox_orf_z, color = 'black')
            
            # plotting the fu orifices
            
            for k in range(0, num_fu_per):
                psi = ((2 * np.pi) / num_fu_per) * k
                fu_x = (cup_ring_dia / 2) * np.cos(i * (2 * np.pi / (num_ox_core - 1))) + ((cup_inset / (np.tan(fu_angle) * 2)) + (cup_dia / 2)) * np.cos(k * (2 * np.pi / num_fu_per))
                fu_y = (cup_ring_dia / 2) * np.sin(i * (2 * np.pi / (num_ox_core - 1))) + ((cup_inset / (np.tan(fu_angle) * 2)) + (cup_dia / 2)) * np.sin(k * (2 * np.pi / num_fu_per))
                center_point = np.array([fu_x, fu_y, -(cup_inset / 2)])
                direction_vector = np.array([np.cos(psi) * np.cos(fu_angle), np.sin(psi) * np.cos(fu_angle), np.sin(fu_angle),])
                
                self.plot_rotating_circle_3d(fu_dia / 2, direction_vector, center_point)
                
        for k in range(0, num_fu_center):
            psi = ((2 * np.pi) / num_fu_center) * k
            center_center_pt = np.array([((cup_inset / (np.tan(fu_angle) * 2)) + (cup_dia / 2)) * np.cos(k * (2 * np.pi / num_fu_center)), ((cup_inset / (np.tan(fu_angle) * 2)) + (cup_dia / 2)) * np.sin(k * (2 * np.pi / num_fu_center)), -(cup_inset / 2)])
            center_direction_vector = np.array([np.cos(psi) * np.cos(fu_angle), np.sin(psi) * np.cos(fu_angle), np.sin(fu_angle)])
            self.plot_rotating_circle_3d(fu_dia / 2, center_direction_vector, center_center_pt)
        
    
    
        # plotting the film cooling orifices
        
        for k in range(num_film):
            film_ring_rad = film_ring_dia / 2
            offset_x = film_ring_rad * np.cos(k * (2 * np.pi) / (num_film))
            offset_y = film_ring_rad * np.sin(k * (2 * np.pi) / (num_film))
            
            film_xs = (film_dia / 2) * np.cos(theta) + offset_x
            film_ys = (film_dia / 2) * np.sin(theta) + offset_y
            film_zs = 0 # flush with the injector face
            
            plt.plot(film_xs, film_ys, film_zs, color = 'black')
    
        
    def go(self):
        
        t_incr = 1 / self.res
        time = 0
        
        while(time < .01):
            # physics section of the code
            for fu in self.fu_list:
                fu.pos_vectors.append(fu.vel_vectors[-1] * t_incr + fu.pos_vectors[-1])
                fu.vel_vectors.append(fu.acc_vectors[-1] * t_incr + fu.vel_vectors[-1])
                fu.acc_vectors.append(np.array([0, 0, 9.8]))


            for ox in self.ox_list:
                
                ox.pos_vectors.append(ox.vel_vectors[-1] * t_incr + ox.pos_vectors[-1])
                ox.vel_vectors.append(ox.acc_vectors[-1] * t_incr + ox.vel_vectors[-1])
                ox.acc_vectors.append(np.array([0, 0, 9.8]))

            
            
            
            
            time += t_incr

        
    #TODO run the simulation without plotting


    def disp(self):
        
        for fu in self.fu_list:
            fu.plot_simple()
        
        for ox in self.ox_list:
            ox.plot_simple('red')
        
        plt.show()
    #TODO plot the results from the simulation on the 3D graph
