import numpy as np
# import matplotlib.pyplot as plt
from plot import plot_vector, plot_vertice


class VerticeAnalysis:

    def __init__(self, p_x, p_y, mag_x, mag_y):
        self.v_x = None
        self.v_y = None
        self.vertice = None
        self.vertice_q = None
        self.charge = None
        self.x = p_x
        self.y = p_y
        self.mx = mag_x
        self.my = mag_y
        self.Ns = len(self.x)
        self.Spin = np.random.choice([-1, 1], size=self.Ns)
        self._Lx = max(self.x) - min(self.x)
        self._Ly = max(self.y) - min(self.y) + np.sin(np.deg2rad(60))
        self._calculate_vertices_()
        self._calculate_vertices_q_()
        self._get_charge()

    def _calculate_vertices_(self):
        """
        Calculate the vertices of the model.
        """
        vx0_k = [0.5, 0.5]
        vy0_k = [0.5 * np.tan(np.deg2rad(30)), np.tan(np.deg2rad(60.)) - 0.5 * np.tan(np.deg2rad(30))]
        vx0_t = 1.5
        vy0_t = 0.5 * np.tan(np.deg2rad(60.))
        v_x = []
        v_y = []
        nx = int(np.sqrt(self.Ns / 3))
        ny = int(np.sqrt(self.Ns / 3))

        for j in range(ny):
            for i in range(nx):
                for k in range(2):
                    v_x.append(vx0_k[k] + 2 * i + j % 2)
                    v_y.append(vy0_k[k] + 2 * np.sin(np.deg2rad(60.)) * j)
                v_x.append(vx0_t + 2 * i + j % 2)
                v_y.append(vy0_t + 2 * np.sin(np.deg2rad(60)) * j)

        v_x = np.array(v_x) - 0.5 * self._Lx
        v_y = np.array(v_y) - 0.5 * self._Ly
        self.v_x = v_x
        self.v_y = v_y

    def _calculate_vertices_q_(self):
        """
        Calculates the vertices charge of the model.
        """
        vertice = []
        vertice_q = []

        for i in range(len(self.v_x)):
            ver = []
            ver_q = []
            for ni in np.arange(-2, 3):
                for nj in np.arange(-2, 3):
                    for j in range(len(self.x)):
                        dx = self.v_x[i] - self.x[j] + float(ni * self._Lx)
                        dy = self.v_y[i] - self.y[j] + float(nj * self._Ly)
                        dist = np.sqrt(dx ** 2 + dy ** 2)
                        if (i + 1) % 3 != 0:
                            '''
                            Kagome
                            '''
                            if 1.1 > dist > 0.1:
                                lij: float = dx * self.mx[j] + dy * self.my[j]
                                ver.append(j)
                                ver_q.append(lij)
                        else:
                            '''
                            Triangular
                            '''
                            if 1.3 > dist > 0.9:
                                lij = dx * self.mx[j] + dy * self.my[j]
                                ver.append(j)
                                ver_q.append(lij)
            if len(ver) != 0:
                vertice.append(ver)
                vertice_q.append(ver_q)

        self.vertice = vertice
        self.vertice_q = vertice_q

    def _get_charge(self):
        self.charge = np.zeros(len(self.vertice))

        for i in range(len(self.vertice)):
            Q = 0
            k = 0
            for j in self.vertice[i]:
                Q += self.vertice_q[i][k] * self.Spin[j]
                k += 1

            self.charge[i] = Q

    def plot_model(self, **kwargs):
        pos = np.column_stack([self.x, self.y])
        vec = np.column_stack([self.mx * self.Spin, self.my * self.Spin])
        return plot_vector(pos, vec, **kwargs)

    def plot_model_vertice(self, label=False, ax=None):
        self._get_charge()
        return plot_vertice(self.v_x, self.v_y, self.charge, ax=ax, label=label)
