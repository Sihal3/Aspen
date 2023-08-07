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

    def plot_complex(self, res, color = 'red'):
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
        
        for i in range(10, len(self.pos_vectors), res):
            # plot_rotating_circle_3d_2(self.dia[i]/2, self.vel_vectors[i], self.pos_vectors[i], color)
            pass