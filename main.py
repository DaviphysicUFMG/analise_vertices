import numpy as np
# import matplotlib.pyplot as plt
from plot import plot_vector, plot_vertice


class VerticeAnalysis:

    def __init__(self, p_x, p_y, mag_x, mag_y, spin=None, nx=50, ny=50):
        self.v_x = None
        self.v_y = None
        self.vertex = None
        self.vertex_q = None
        self.charge = None
        self.x = p_x
        self.y = p_y
        self.mx = mag_x
        self.my = mag_y
        self.nx = nx
        self.ny = ny
        self.Bx = np.zeros(shape=(nx, ny))
        self.By = np.zeros(shape=(nx, ny))
        self.XX = np.zeros(shape=(nx, ny))
        self.YY = np.zeros(shape=(nx, ny))
        self.Ns = len(self.x)
        if spin is None:
            self.Spin = np.ones(self.Ns)  # np.random.choice([-1, 1], size=self.Ns)
        else:
            self.Spin = spin
        self._Lx = max(self.x) - min(self.x)
        self._Ly = max(self.y) - min(self.y) + np.sin(np.deg2rad(60))
        self._calculate_vertices_()
        self._calculate_vertices_q_()
        self.get_charge()
        self._init_grid()
        self.calc_field()

    def _calculate_vertices_(self):
        """
        Calculate the vertices of the model.
        """
        xmin = min(self.x)
        ymin = min(self.y)
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
                    v_x.append(xmin + vx0_k[k] + 2 * i + j % 2)
                    v_y.append(ymin + vy0_k[k] + 2 * np.sin(np.deg2rad(60.)) * j)
                v_x.append(xmin + vx0_t + 2 * i + j % 2)
                v_y.append(ymin + vy0_t + 2 * np.sin(np.deg2rad(60)) * j)

        v_x = np.array(v_x)
        v_y = np.array(v_y)
        # v_x += 0.5 * self._Lx
        # v_y += 0.5 * self._Ly
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
            for ni in np.arange(-1, 2):
                for nj in np.arange(-1, 2   ):
                    for j in range(len(self.x)):
                        dx = self.v_x[i] - self.x[j] + float(ni * self._Lx)
                        dy = self.v_y[i] - self.y[j] + float(nj * self._Ly)
                        dist = np.sqrt(dx ** 2 + dy ** 2)
                        if (i + 1) % 3 != 0:
                            '''
                            Kagome
                            '''
                            if 1.1 > dist > 0.1:
                                lij = (dx * self.mx[j] + dy * self.my[j]) / dist ** 3
                                ver.append(j)
                                ver_q.append(lij)
                        else:
                            '''
                            Triangular
                            '''
                            if 1.3 > dist > 0.9:
                                lij = (dx * self.mx[j] + dy * self.my[j]) / dist ** 3
                                ver.append(j)
                                ver_q.append(lij)
            if len(ver) != 0:
                vertice.append(ver)
                vertice_q.append(ver_q)

        self.vertex = vertice
        self.vertex_q = vertice_q

    def _init_grid(self):
        xmin, xmax = min(self.x), max(self.x)
        ymin, ymax = min(self.y), max(self.y)
        X = np.linspace(xmin, xmax, self.nx)
        Y = np.linspace(ymin, ymax, self.ny)
        self.XX, self.YY = np.meshgrid(X, Y)

    def calc_field(self):
        self.Bx = np.zeros(shape=(self.nx, self.ny))
        self.By = np.zeros(shape=(self.nx, self.ny))

        for jB in range(self.ny):
            for iB in range(self.nx):
                for i in range(self.Ns):
                    dx = self.XX[iB, jB] - self.x[i]
                    dy = self.YY[iB, jB] - self.y[i]
                    dist = np.sqrt(dx ** 2 + dy ** 2)
                    if dist > 0.1:
                        dx /= dist
                        dy /= dist
                        a1 = 3*((self.mx[i] * dx) + (self.my[i] * dy))
                        self.Bx[iB, jB] += self.Spin[i] * (a1 * dx - self.mx[i]) / dist ** 3
                        self.By[iB, jB] += self.Spin[i] * (a1 * dy - self.my[i]) / dist ** 3

        self.Bx /= self.Ns
        self.By /= self.Ns

    def get_charge(self):
        self.charge = np.zeros(len(self.vertex))

        for i in range(len(self.vertex)):
            Q = 0
            k = 0
            for j in self.vertex[i]:
                Q += self.vertex_q[i][k] * self.Spin[j]
                k += 1

            self.charge[i] = Q

    def plot_model(self, **kwargs):
        pos = np.column_stack([self.x, self.y])
        vec = np.column_stack([self.mx * self.Spin, self.my * self.Spin])
        return plot_vector(pos, vec, **kwargs)

    def plot_model_vertice(self, label=False, ax=None):
        self.get_charge()
        return plot_vertice(self.v_x, self.v_y, self.charge, ax=ax, label=label)
