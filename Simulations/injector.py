import math
from pint import UnitRegistry
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt



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

massflow = 1 * ureg('kilogram / second')
B = 1.2
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

def within_crossflow(fu_jet, list_ox_jets):
    ans = False
    for ox in list_ox_jets:
        x_f = fu_jet.pos_vectors[-1][0]
        x_o = ox.pos_vectors[-1][0]
        y_f = fu_jet.pos_vectors[-1][1]
        y_o = ox.pos_vectors[-1][1]
        z_f = fu_jet.pos_vectors[-1][2]
        z_o = ox.pos_vectors[-1][2]
      
        if (ox.dia[-1] / 2) >= np.abs( np.sqrt ( (x_f - x_o)**2 + (y_f - y_o)**2)):
            ans = ox
            break
    return ans
        
        
class Jet:
    
    
    def __init__(self, num_ox_core, num_fu_per, cup_inset, cup_ring_dia, cup_dia):
        self.num_ox_core = num_ox_core
        self.num_fu_per = num_fu_per
        self.cup_inset = cup_inset
        self.cup_ring_dia = cup_ring_dia
        self.cup_dia = cup_dia
        
        self.pos_vectors = []
        self.dia = []
        self.vel_vectors = []
        self.ab = []

    def plot_simple(self, color= 'blue'):
        xs = []
        ys = []
        zs = []
        
        for pos in self.pos_vectors:
            xs.append(pos[0])
        for pos in self.pos_vectors:
            ys.append(pos[1])
        for pos in self.pos_vectors:
            zs.append(pos[2])

        plt.plot(xs, ys, zs, color = color)


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
        return v
    return v / norm

class Ox_Jet(Jet):
    def __init__(self, id_num, num_ox_core, num_fu_per, cup_inset, cup_ring_dia, cup_dia, ox_dia):
        super().__init__(num_ox_core, num_fu_per, cup_inset, cup_ring_dia, cup_dia)
        self.id_num = id_num
        self.ox_dia = ox_dia
        self.dia.append(self.ox_dia)
        
    
    def set_origins(self):
        # sets the initial positions for the jets 
        cup_ring_rad = self.cup_ring_dia / 2
        if self.id_num == 0:
            pos_x = 0
            pos_y = 0
        else:
            pos_x = cup_ring_rad * np.cos(self.id_num * (2 * np.pi) / (self.num_ox_core - 1))
            pos_y = cup_ring_rad * np.sin(self.id_num * (2 * np.pi) / (self.num_ox_core - 1))

        pos_z = - self.cup_inset
            
        self.pos_vectors.append(np.array([pos_x, pos_y, pos_z]))
    
    def set_initial_velocity(self, vel):
        # takes in a scalar velocity and makes it into a vector assuming it is parallel to chamber
        self.vel_vectors.append(np.array([0, 0, vel]))
    
        
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
            ox_location = math.floor(self.id_num / self.num_fu_per)
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
        print(f"phi {np.degrees(psi)} fu {np.degrees(ang)}")
        direction_vector = np.array([np.cos(psi) * np.cos(ang), np.sin(psi) * np.cos(ang), np.sin(ang)])
        unit_vel = vel * normalize(direction_vector)
        self.vel_vectors.append(unit_vel)
        
        
        
    
class Sim:
    def __init__(self):
        self.fu_list = []
        self.ox_list = []
        
        self.chamber_dia = -1
        self.chamber_len = -1
        self.cup_ring_dia = -1
        self.cup_dia = -1
        self.fu_inj_angle = -1
        self.fu_inj_dia = -1
        self.ox_inj_dia = -1
        self.cup_inset = -1
        self.num_ox_core = -1
        self.num_fu_per = -1
        self.Cd = .8
        
    def area_ellipse(self, a, b):
        return a * b * np.pi
    
    def y_vel_crossflow(self, t, fu_jet, ox_jet):
        # adds the new velocity for the fuel jet
        # also updates the direction for the jet
        v_x = fu_jet.vel_vectors[-1][0]
        v_y = fu_jet.vel_vectors[-1][1]
        delta_v_z = ox_jet.vel_vectors[-1][2] - fu_jet.vel_vectors[-1][2]
        delta_v_z *= ureg('in/s')
        delta_v_z.ito('m/s')
        num = 2 * np.cos(self.fu_inj_angle) * self.ox_density * delta_v_z ** 2 * self.Cd * t
        dia = fu_jet.dia[-1]  * ureg('in')
        dia.ito('m')
        den = self.fu_density * np.pi * dia
        resultant = ((num / den).magnitude) * ureg('m/s')
        resultant.ito('in/s')
        fu_jet.vel_vectors.append(np.array([v_x, v_y, (resultant.magnitude) + fu_jet.vel_vectors[-1][2]]))
    
    def chamber_plot(self, chamber_dia_u, chamber_len_u, num_ox_core, cup_ring_dia_u, num_fu_per, cup_dia_u, fu_inj_angle_u, fu_inj_dia_u, ox_inj_dia_u, cup_inset_u):
    
        # strip all of the inputted dimensions of their units
        chamber_dia = chamber_dia_u.magnitude
        chamber_len = chamber_len_u.magnitude
        cup_ring_dia = cup_ring_dia_u.magnitude
        cup_dia = cup_dia_u.magnitude
        fu_inj_angle = fu_inj_angle_u.magnitude
        fu_inj_dia = fu_inj_dia_u.magnitude
        ox_inj_dia = ox_inj_dia_u.magnitude
        cup_inset = cup_inset_u.magnitude
        
        fu_inj_angle = np.radians(fu_inj_angle)
        
        ax.set_xlim(-10, 10) 
        ax.set_ylim(-10, 10)
        ax.set_zlim(-10, 10)
        
        r = chamber_dia / 2
        
        theta = np.linspace(0, 2 * np.pi, 100)
        injector_face_xs = r * np.cos(theta)
        injector_face_ys = r * np.sin(theta)
        
        
        
        plt.plot(injector_face_xs, injector_face_ys, 0, color = 'black')
        plt.plot(injector_face_xs, injector_face_ys, chamber_len, color = 'black')

        cup_z = -cup_inset 
        
        # center injector
            # outline flush with the injector
        center_cup_flush_dia = cup_dia + 2 * (cup_inset / np.tan(fu_inj_angle))
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
        center_ox_orf_xs = (ox_inj_dia / 2) * np.cos(theta)
        center_ox_orf_ys = (ox_inj_dia / 2) * np.sin(theta)
        center_ox_orf_z = - cup_inset
        
        plt.plot(center_ox_orf_xs, center_ox_orf_ys, center_ox_orf_z, color = 'black')
        
        
        for i in range(0, num_ox_core - 1):
            # plot the 3 circles representing the cups 
            
            # the offset for the cups that creates the cup pattern
            cup_ring_rad = cup_ring_dia / 2
            offset_x = cup_ring_rad * np.cos(i * (2 * np.pi) / (num_ox_core - 1))
            offset_y = cup_ring_rad * np.sin(i * (2 * np.pi) / (num_ox_core - 1))
            
            # the first one is the outline flush with the injector face
            
            cup_flush_dia = cup_dia + 2 * (cup_inset / np.tan(fu_inj_angle))
            
            
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
            ox_orf_xs = (ox_inj_dia / 2) * np.cos(theta)
            ox_orf_ys = (ox_inj_dia / 2) * np.sin(theta)
            ox_orf_z = - cup_inset
            
            ox_orf_xs += offset_x
            ox_orf_ys += offset_y
            
            plt.plot(ox_orf_xs, ox_orf_ys, ox_orf_z, color = 'black')
            
            
            # plotting of the fuel jet orifices
            
            
            
            for k in range(0, num_fu_per):
                psi = ((2 * np.pi) / num_fu_per) * k
                fu_x = (cup_ring_dia / 2) * np.cos(i * (2 * np.pi / (num_ox_core - 1))) + ((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.cos(k * (2 * np.pi / num_fu_per))
                fu_y = (cup_ring_dia / 2) * np.sin(i * (2 * np.pi / (num_ox_core - 1))) + ((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.sin(k * (2 * np.pi / num_fu_per))
                center_point = np.array([fu_x, fu_y, -(cup_inset / 2)])
                direction_vector = np.array([np.cos(psi) * np.cos(fu_inj_angle), np.sin(psi) * np.cos(fu_inj_angle), np.sin(fu_inj_angle),])
            
                plot_rotating_circle_3d(fu_inj_dia / 2, direction_vector, center_point)
                
        for k in range(0, num_fu_per):
            psi = ((2 * np.pi) / num_fu_per) * k
            center_center_pt = np.array([((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.cos(k * (2 * np.pi / num_fu_per)), ((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.sin(k * (2 * np.pi / num_fu_per)), -(cup_inset / 2)])
            center_direction_vector = np.array([np.cos(psi) * np.cos(fu_inj_angle), np.sin(psi) * np.cos(fu_inj_angle), np.sin(fu_inj_angle)])
            plot_rotating_circle_3d(fu_inj_dia / 2, center_direction_vector, center_center_pt)
    
    def create_scenario(self, chamber_dia_u, chamber_len_u, num_ox_core, cup_ring_dia_u, num_fu_per, cup_dia_u, fu_inj_angle_u, fu_inj_dia_u, ox_inj_dia_u, cup_inset_u, fu_vel_u, ox_vel_u, ox_density, fu_density):
        # strip all of the inputted dimensions of their units
        self.chamber_dia = chamber_dia_u.magnitude
        self.chamber_len = chamber_len_u.magnitude
        self.cup_ring_dia = cup_ring_dia_u.magnitude
        self.cup_dia = cup_dia_u.magnitude
        fu_inj_angle = fu_inj_angle_u.magnitude
        self.fu_inj_dia = fu_inj_dia_u.magnitude
        self.ox_inj_dia = ox_inj_dia_u.magnitude
        self.cup_inset = cup_inset_u.magnitude
        self.num_ox_core = num_ox_core
        self.num_fu_per = num_fu_per
        self.fu_density = fu_density
        self.ox_density = ox_density
        

        self.fu_inj_angle = np.radians(fu_inj_angle)
        fu_vel = fu_vel_u.magnitude
        ox_vel = ox_vel_u.magnitude
        
        
        
        num_fu = num_ox_core * num_fu_per
        num_ox = num_ox_core
        
        
        #! Fill the lists of jets
        for fu_id in range(num_fu):
            # create Fu_Jet obj and add to list 
            new_fu_jet = Fu_Jet(fu_id, self.num_ox_core, self.num_fu_per, self.cup_inset, self.cup_ring_dia, self.cup_dia, self.fu_inj_dia, self.fu_inj_angle) 
            new_fu_jet.set_initial_velocity(fu_vel)
            new_fu_jet.set_origins()
            self.fu_list.append(new_fu_jet)
        
        for ox_id in range(num_ox):
            # create Ox_Jet obj and add to list
            new_ox_jet = Ox_Jet(ox_id, self.num_ox_core, self.num_fu_per, self.cup_inset, self.cup_ring_dia, self.cup_dia, self.ox_inj_dia)
            new_ox_jet.set_origins()
            new_ox_jet.set_initial_velocity(ox_vel)
            self.ox_list.append(new_ox_jet)
        
        print(self.fu_list[2].pos_vectors)
    
    def go():
        pass 
    
    def simple(self):
        # Simplest model, just plot straight lines from each of the orifices
        time = 0
        t_incr = 0.00002 
        
        while time < .1:
            # start with the ox jets
            for ox in self.ox_list:
                ox.pos_vectors.append(ox.vel_vectors[-1] * t_incr + ox.pos_vectors[-1])
                ox.vel_vectors.append(np.array([0,0,0]) * t_incr + ox.vel_vectors[-1]) # no acceleration or change in trajectory
            
            # then work with the fu jets 
            for fu in self.fu_list:
                fu.pos_vectors.append(fu.vel_vectors[-1] * t_incr + fu.pos_vectors[-1])
                
                if within_crossflow(fu, self.ox_list):
                    
                    ox_in_question = within_crossflow(fu, self.ox_list)
                    self.y_vel_crossflow(t_incr, fu, ox_in_question)
                
            time += t_incr
        
        for ox in self.ox_list:
            ox.plot_simple('red')
        for fu in self.fu_list:
            fu.plot_simple()

    def show():
        pass

def cross_flow_region(Cd, p_g, p_l, u_l, u_g, d_fu, fu_angle, x):
    fu_angle = math.radians(fu_angle)
    u_y = u_l * math.sin(fu_angle)
    u_x = u_l * math.cos(fu_angle)
    
    delta_u = u_g - u_y
    
    num = Cd * p_g * delta_u ** 2 * x ** 2 * math.cos(fu_angle) 
    den = math.pi * p_l * u_x ** 2 * d_fu
    
    vel = ( u_y * x ) / u_x
    
    return num / den + vel


def pre_cross_flow_region(u_l, fu_angle, x):
    fu_angle = math.radians(fu_angle)
    u_y = u_l * math.sin(fu_angle)
    u_x = u_l * math.cos(fu_angle)
    
    return (u_y /u_x) * x

def post_cross_flow_region(Cd, p_g, p_l, u_l, u_g, d_fu, fu_angle, x, u_y0):
    fu_angle = math.radians(fu_angle)
    u_y = u_l * math.sin(fu_angle)
    u_x = u_l * math.cos(fu_angle)
    delta_u = u_g - u_y
    
    num = 2 * math.cos(fu_angle) * p_g * delta_u ** 2 * Cd * x
    den = p_l * math.pi * d_fu * u_x
    print(( (num / den + u_y0) / u_x ))
    return  ( (num / den + u_y0) / u_x ) * x


def x_1_locator(cup_dia, d_ox, cup_inset, fu_angle):
    fu_angle = math.radians(fu_angle)
    
    return (cup_dia - d_ox + cup_inset * math.tan(fu_angle)) / 2





def plot_model(Cd, p_g, p_l, u_l, u_g, d_fu, fu_angle, cup_dia, d_ox, cup_inset, dx):
    
    
    x_1 = x_1_locator(cup_dia, d_ox, cup_inset, fu_angle)
    ys = []
    xs = []
    C = 1 / dx

    i = 0
    while i < x_1.magnitude:
        ys.append(pre_cross_flow_region(u_l, fu_angle, i * ureg("m")).magnitude)
        xs.append(i)
        i += dx

    while i < x_1.magnitude + d_ox.magnitude:
        ys.append(cross_flow_region(Cd, p_g, p_l, u_l, u_g, d_fu, fu_angle, i * ureg("m")).magnitude)
        xs.append(i)
        i += dx

    u_y0 = u_l * math.sin(fu_angle)

    while i < 2 * (x_1.magnitude + d_ox.magnitude):
    
        ys.append(post_cross_flow_region(Cd, p_g, p_l, u_l, u_g, d_fu, fu_angle, i * ureg("m"), u_y0).magnitude)
        xs.append(i)
        i += dx

    plt.plot(xs, ys)
    plt.axvline(x = x_1.magnitude, color = 'b', label = 'axvline - full height')
    plt.axvline(x = (x_1 + d_ox).magnitude, color = 'r', label = 'axvline - full height')
    plt.show()


plt.rcParams['figure.figsize'] = ( 10,10)
ax = plt.axes(projection='3d')



def plot_rotating_ellipse_3d(a, b, direction_vector, center_point, theta):
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
    ax.plot(x, y, z, color='b', alpha=0.6)


def plot_rotating_circle_3d(radius, direction_vector, center):
    # radius is int
    # direction_vector is an np array
    # center is an np array representing a point
    
    # Normalize the direction vector to ensure it is a unit vector
    direction_vector = direction_vector / np.linalg.norm(direction_vector)
    
    # Generate points on the circle in 3D using parametric equations
    t = np.linspace(0, 2 * np.pi, 100)
    # Generate the normal vector to the plane of the circle
    n = np.cross(direction_vector, np.array([0, 0, 1]))
    n = n / np.linalg.norm(n)
    
    # Generate the circle points in 3D space using the direction vector and normal vector
    x = center[0] + radius * (np.cos(t) * direction_vector[0] + np.sin(t) * n[0])
    y = center[1] + radius * (np.cos(t) * direction_vector[1] + np.sin(t) * n[1])
    z = center[2] + radius * (np.cos(t) * direction_vector[2] + np.sin(t) * n[2])
    
    
    
    # Plot the circle
    ax.plot(x, y, z, color='black', alpha=0.6)





def chamber_plot(chamber_dia_u, chamber_len_u, num_ox_core, cup_ring_dia_u, num_fu_per, cup_dia_u, fu_inj_angle_u, fu_inj_dia_u, ox_inj_dia_u, cup_inset_u):
    
    # strip all of the inputted dimensions of their units
    chamber_dia = chamber_dia_u.magnitude
    chamber_len = chamber_len_u.magnitude
    cup_ring_dia = cup_ring_dia_u.magnitude
    cup_dia = cup_dia_u.magnitude
    fu_inj_angle = fu_inj_angle_u.magnitude
    fu_inj_dia = fu_inj_dia_u.magnitude
    ox_inj_dia = ox_inj_dia_u.magnitude
    cup_inset = cup_inset_u.magnitude
    
    fu_inj_angle = np.radians(fu_inj_angle)
    
    ax.set_xlim(-10, 10) 
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    
    r = chamber_dia / 2
    
    theta = np.linspace(0, 2 * np.pi, 100)
    injector_face_xs = r * np.cos(theta)
    injector_face_ys = r * np.sin(theta)
    
    
    
    plt.plot(injector_face_xs, injector_face_ys, 0, color = 'black')
    plt.plot(injector_face_xs, injector_face_ys, chamber_len, color = 'black')

    cup_z = -cup_inset 
    
    # center injector
        # outline flush with the injector
    center_cup_flush_dia = cup_dia + 2 * (cup_inset / np.tan(fu_inj_angle))
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
    center_ox_orf_xs = (ox_inj_dia / 2) * np.cos(theta)
    center_ox_orf_ys = (ox_inj_dia / 2) * np.sin(theta)
    center_ox_orf_z = - cup_inset
    
    plt.plot(center_ox_orf_xs, center_ox_orf_ys, center_ox_orf_z, color = 'black')
    
    
    for i in range(0, num_ox_core - 1):
        # plot the 3 circles representing the cups 
        
        # the offset for the cups that creates the cup pattern
        cup_ring_rad = cup_ring_dia / 2
        offset_x = cup_ring_rad * np.cos(i * (2 * np.pi) / (num_ox_core - 1))
        offset_y = cup_ring_rad * np.sin(i * (2 * np.pi) / (num_ox_core - 1))
        
        # the first one is the outline flush with the injector face
        
        cup_flush_dia = cup_dia + 2 * (cup_inset / np.tan(fu_inj_angle))
        
        
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
        ox_orf_xs = (ox_inj_dia / 2) * np.cos(theta)
        ox_orf_ys = (ox_inj_dia / 2) * np.sin(theta)
        ox_orf_z = - cup_inset
        
        ox_orf_xs += offset_x
        ox_orf_ys += offset_y
        
        plt.plot(ox_orf_xs, ox_orf_ys, ox_orf_z, color = 'black')
        
        
        # plotting of the fuel jet orifices
        
        
        
        for k in range(0, num_fu_per):
            psi = ((2 * np.pi) / num_fu_per) * k
            fu_x = (cup_ring_dia / 2) * np.cos(i * (2 * np.pi / (num_ox_core - 1))) + ((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.cos(k * (2 * np.pi / num_fu_per))
            fu_y = (cup_ring_dia / 2) * np.sin(i * (2 * np.pi / (num_ox_core - 1))) + ((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.sin(k * (2 * np.pi / num_fu_per))
            center_point = np.array([fu_x, fu_y, -(cup_inset / 2)])
            direction_vector = np.array([np.cos(psi) * np.cos(fu_inj_angle), np.sin(psi) * np.cos(fu_inj_angle), np.sin(fu_inj_angle),])
        
            plot_rotating_circle_3d(fu_inj_dia / 2, direction_vector, center_point)
            
    for k in range(0, num_fu_per):
        psi = ((2 * np.pi) / num_fu_per) * k
        center_center_pt = np.array([((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.cos(k * (2 * np.pi / num_fu_per)), ((cup_inset / (np.tan(fu_inj_angle) * 2)) + (cup_dia / 2)) * np.sin(k * (2 * np.pi / num_fu_per)), -(cup_inset / 2)])
        center_direction_vector = np.array([np.cos(psi) * np.cos(fu_inj_angle), np.sin(psi) * np.cos(fu_inj_angle), np.sin(fu_inj_angle)])
        plot_rotating_circle_3d(fu_inj_dia / 2, center_direction_vector, center_center_pt)


def plot_ox(num_ox_core, cup_ring_dia_u, cup_inset_u, chamber_len_u):
    cup_ring_dia = cup_ring_dia_u.magnitude
    cup_inset = cup_inset_u.magnitude
    chamber_len = chamber_len_u.magnitude
    # first step is getting an array of the x y z positions of the middle of each gas jet
    
    centers = []
    
    center_z = - cup_inset
    for i in range(0, num_ox_core - 1):
        # plot the 3 circles representing the cups 
        
        # the offset for the cups that creates the cup pattern
        cup_ring_rad = cup_ring_dia / 2
        center_x = cup_ring_rad * np.cos(i * (2 * np.pi) / (num_ox_core - 1))
        center_y = cup_ring_rad * np.sin(i * (2 * np.pi) / (num_ox_core - 1))
        
        center_position = [center_x, center_y, center_z]
        centers.append(center_position)

    
    # test
    for i in centers:
        
        xs = np.ones(10) * i[0]
        
        ys = np.ones(10) * i[1]
        zs = i[2] + 1 * np.linspace(0,10,10)
        
        xs = np.array(xs)
        ys = np.array(ys)
        zs = np.array(zs)
        print(xs , ys, zs)
        plt.plot(xs,ys,zs, color='red')
        



def inj_func(mDot, B, chamber_pressure, ign_chamber_press):
    print()
    drill_dia = {'1/64' : 0.396875 * ureg('mm'), '1/32' : 0.79375 * ureg('mm'), '1/16' : 1.5875 * ureg('mm')}
    
    ox_cores = 7
    fu_inj_per = 7
    film_cooling_orifices = 12
    fu_inj_angle = 50
    
    fu_temp = 290 * ureg('kelvin')
    ox_temp = 250.15 * ureg('kelvin')

    fu_v = 20 * ureg('in/sec')
    ox_v = 120 * ureg('in/sec')
    
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
    
    chamber_plot(3.5 * ureg('inch'), 12 * ureg('inch'), ox_cores, 2 * ureg('inch'), fu_inj_per, element_ox_dia.to(ureg.inch) * 1.5, fu_inj_angle * ureg('degree'), element_fu_dia.to(ureg.inch), element_ox_dia.to(ureg.inch), .25 * ureg('inch'))
    
    test = Sim()
    test.create_scenario(3.5 * ureg('inch'), 12 * ureg('inch'), ox_cores, 2 * ureg('inch'), fu_inj_per, element_ox_dia.to(ureg.inch) * 1.5, fu_inj_angle * ureg('degree'), element_fu_dia.to(ureg.inch), element_ox_dia.to(ureg.inch), .25 * ureg('inch'), fu_v, ox_v, ox_density, fu_density)
    test.simple()
    


inj_func(massflow, B, chamber_pressure, chamber_pressure * 1.2)

# chamber_plot(3.6, 12, 7, 2.2, 7, .5, np.pi /4 , .05, .25, .2)
plt.show()