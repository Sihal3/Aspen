from jet import Jet
import numpy as np
import CoolProp.CoolProp as CP


class Fu_Jet(Jet):
    def __init__(self, id_num, package):
        super().__init__(package)
        self.package = package
        self.fu_angle_u = package['fu_angle_u']
        self.id_num = id_num
        self.dia.append(package['fu_dia_u'].magnitude)
        self.fu_temp = package['fu_temp_u'].magnitude
        self.fuel = package['fuel']
        self.inj_pressure = package['inj_pressure_u'].magnitude
        self.atomized = False
        self.cbl = 0 # add to this as it passes through the ox jet
        
    def add_phase(self):
        phase = CP.PhaseSI("P", self.inj_pressure, "T", self.fu_temp, self.fuel)
        
        if phase == "liquid" and not self.atomized:
            self.phase.append(1)
        else:
            self.phase.append(0) # assuming it can either be gas or liq

    def atomize_it(self):
        self.atomized = True
    
    
    def add_dia(self, dia):
        if self.atomized:
            self.dia.append(self.dia[-1] * 1.01)
        else:
            super().add_dia(dia)
        
    

    def add_density(self):
        den = CP.PropsSI('D', 'T', self.temp[-1], 'P', self.package['chamber_pressure_u'].magnitude, self.fuel)
        self.density.append(den)

    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm

    def set_origins(self):
        
        if self.id_num < self.num_fu_center:

            pos_x = ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.cos(self.id_num * (2 * np.pi / self.num_fu_center))
            pos_y = ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.sin(self.id_num * (2 * np.pi / self.num_fu_center))

        
        else:
            ox_location = np.floor(self.id_num / self.num_fu_per)
            fu_location = self.id_num % self.num_fu_per
            
            pos_x = (self.cup_ring_dia_u / 2) * np.cos(ox_location * (2 * np.pi / (self.num_ox_core - 1))) + ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.cos(fu_location * (2 * np.pi / self.num_fu_per))
            pos_y = (self.cup_ring_dia_u / 2) * np.sin(ox_location * (2 * np.pi / (self.num_ox_core - 1))) + ((self.cup_inset_u / (np.tan(self.fu_angle_u) * 2)) + (self.cup_dia_u / 2)) * np.sin(fu_location * (2 * np.pi / self.num_fu_per))
                
        pos_z = -(self.cup_inset_u / 2)        
                
        self.pos_vectors.append(np.array([pos_x.magnitude, pos_y.magnitude, pos_z.magnitude]))

    def compound_cbl(self, t_incr):
        v_x = self.vel_vectors[-1][0]
        v_y = self.vel_vectors[-1][1]
        p_x = v_x * t_incr
        p_y = v_y * t_incr
        
        dist_in = np.sqrt((p_x ** 2) + (p_y ** 2))
        self.cbl += dist_in
        
    
    def set_initial_velocity(self):
        # set vector for velocity
        if self.id_num < self.num_fu_center:
            fu_location = self.id_num % self.num_fu_center
            psi = (((2 * np.pi) / self.num_fu_center) * fu_location)
        else:
            fu_location = self.id_num % self.num_fu_per
            psi = (((2 * np.pi) / self.num_fu_per) * fu_location)

        vel = self.package['fu_vel_u'].magnitude
        ang = (self.fu_angle_u - np.radians(270)).magnitude
        direction_vector = np.array([np.cos(psi) * np.cos(ang), np.sin(psi) * np.cos(ang), np.sin(ang)])
        unit_vel = vel * self.normalize(direction_vector)
        self.vel_vectors.append(unit_vel)
        

    def calc_cross_traj(self, ox_jet):
        # calculates the new Z acceleration for the fu jet
        
        # determines if the fu jet is within the crossflow region
        
        theta = (np.pi / 2) - self.fu_angle_u.magnitude
        Cd = 1.67 # taken from 'Atomizations and Sprays'
        num = ox_jet.density[-1] * (ox_jet.vel_vectors[-1][2] - self.vel_vectors[-1][2]) ** 2 * np.cos(theta) * Cd
        den = self.dia[-1] * np.pi * self.density[-1]
        
        accel = num / den
        accel_vector = np.array([0, 0, accel])
        return accel_vector

    def calc_deformation(self, ox_jet):
        #TODO Update the a b values for the jet based on how it is deformed
        pass