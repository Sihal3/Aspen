import matplotlib.pyplot as plt
import numpy as np
from dimensions import ureg


class Jet:
    def __init__(self, package):
        self.num_ox_core = package['num_ox_core']
        self.num_fu_per = package['num_fu_per']
        self.num_fu_center = package['num_fu_center']
        self.cup_inset_u = package['cup_inset_u']
        self.cup_ring_dia_u = package['cup_ring_dia_u']
        self.cup_dia_u = package['cup_dia_u']
        
        self.pos_vectors = []
        self.dia = []
        self.vel_vectors = []
        self.ab = []
        self.acc_vectors = []
        self.phase = [] # 1 = liq, 0 = gas
        self.density = [] 
        self.temp = [] # I know I will eventually want this but thats for later

    

    def plot_rotating_circle_3d_2(self, radius, direction_vector, center_point, index, color = 'blue'):
        
        # Normalize the direction vector to ensure it is a unit vector
        direction_vector = direction_vector / np.linalg.norm(direction_vector)
        
        # Generate points on the circle in 3D using parametric equations
        t = np.linspace(0, 2 * np.pi, 100)
        
        # Calculate a vector orthogonal to the direction vector using cross product
        orthogonal_vector = np.cross(direction_vector, [0, 0, 1])
        
        # If the direction vector is close to (0, 0, 1), use a different orthogonal vector
        if np.allclose(orthogonal_vector, [0, 0, 0]):
            orthogonal_vector = np.cross(direction_vector, [1, 0, 0])
        
        # Normalize the orthogonal vector to ensure it is a unit vector
        orthogonal_vector = orthogonal_vector / np.linalg.norm(orthogonal_vector)
        
        # Calculate the normal vector to the circle plane using cross product
        normal_vector = np.cross(direction_vector, orthogonal_vector)
        
        # Generate the circle points in 3D space using the rotated direction vector and center point
        circle_points = center_point[:, None] + radius * (orthogonal_vector[:, None] * np.cos(t) + normal_vector[:, None] * np.sin(t))
        
        # Create a 3D plot

        # Plot the circle
        line_style = '-'
        if self.phase[index] == -1:
            line_style = '-'
            print("USE SET_PHASE")
        elif self.phase[index] == 1:
            line_style = '-'
        elif self.phase[index] == 0:
            line_style = ":"
        
        plt.plot(circle_points[0], circle_points[1], circle_points[2], linestyle= line_style, color=color, alpha=.6)
        
        # Plot the circle
    def add_temp(self, temp):
        self.temp.append(temp)

    def set_acc(self):
        self.acc_vectors.append(np.array([0, 0, 0]))

    def plot_simple(self, color= 'blue'):
        xs = []
        ys = []
        zs = []
        
        for pos in self.pos_vectors:
            x = pos[0] * ureg('m')
            x.ito(ureg('inch'))
            y = pos[1] * ureg('m')
            y.ito(ureg('inch'))
            z = pos[2] * ureg('m')
            z.ito(ureg('inch'))
            xs.append(x.magnitude)
            ys.append(y.magnitude)
            zs.append(z.magnitude)

        plt.plot(xs, ys, zs, color = color)

    def plot_complex(self, res, color = 'blue'):
        xs = []
        ys = []
        zs = []
        

        
        
        for pos in self.pos_vectors:
            x = pos[0] * ureg('m')
            x.ito(ureg('inch'))
            y = pos[1] * ureg('m')
            y.ito(ureg('inch'))
            z = pos[2] * ureg('m')
            z.ito(ureg('inch'))
            xs.append(x.magnitude)
            ys.append(y.magnitude)
            zs.append(z.magnitude)
        
        
        
        plt.plot(xs, ys, zs, color = color)
        
        for i in range(10, len(self.pos_vectors), res):
            rad = ((self.dia[i] / 2) * ureg('meter')).to('inch').magnitude
            x_in = (self.pos_vectors[i][0] * ureg('meter')).to('inch').magnitude
            y_in = (self.pos_vectors[i][1] * ureg('meter')).to('inch').magnitude
            z_in = (self.pos_vectors[i][2] * ureg('meter')).to('inch').magnitude
            
            pos = np.array([x_in, y_in, z_in])
            
            self.plot_rotating_circle_3d_2(rad, self.vel_vectors[i], pos, i, color)
            
    
    
    def add_dia(self, dia):
        # takes in a diameter as int_u and adds it to the self.dia list
        self.dia.append(dia)