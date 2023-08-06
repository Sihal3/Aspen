import matplotlib.pyplot as plt

class Jet:
    def __init__(self, package):
        self.num_ox_core = package['num_ox_core']
        self.num_fu_per = package['num_fu_per']
        self.num_fu_center = package['num_fu_center']
        self.cup_inset_u = package['cup_inset_u']
        self.cup_ring_dia_u = package['cup_ring_dia_u']
        
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