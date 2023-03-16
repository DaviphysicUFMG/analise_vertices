import numpy as np
from plot import plot_vector, plot_vertice


class VerticeAnalysis:

    def __init__(self, p_x, p_y, mag_x, mag_y, spin=None, nx=50, ny=50):
        """
        Vertice analysis class.

        Parameters
        ----------
        p_x : array_like
            x-coordinates of the vertices.
        p_y : array_like
            y-coordinates of the vertices.
        mag_x : array_like
            x-component of the magnetic field.
        mag_y : array_like
            y-component of the magnetic field.
        spin : array_like, optional
            Spin of the model.
        nx : int, optional
            Number of grid points in the x-direction.
        ny : int, optional
            Number of grid points in the y-direction.
        """
        self.p_x = np.array(p_x, dtype=float)
        self.p_y = np.array(p_y, dtype=float)
        self.mag_x = np.array(mag_x, dtype=float)
        self.mag_y = np.array(mag_y, dtype=float)
        self.spin = None
        self.v_x = None
        self.v_y = None
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
        self._Lx = max(self.x) - min(self.x)
        self._Ly = max(self.y) - min(self.y) + np.sin(np.deg2rad(60))

        self._calc_vertex_()
        self._calc_vertex_charge()

        if spin is None:
            self.Spin = np.ones(self.Ns, dtype=int)  # np.random.choice([-1, 1], size=self.Ns)
        else:
            self.set_spin(spin)

        # self._calculate_vertices_q_()
        # self.get_charge()
        # self._init_grid()
        # self.calc_field()

    def _calc_vertex_(self):
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

        v_x = np.array(v_x, dtype=float)
        v_y = np.array(v_y, dtype=float)
        self.v_x = v_x
        self.v_y = v_y

    def _calc_vertex_charge(self):
        i_vertex = []
        j_vertex = []
        i_vertex_q = []
        for i in range(len(self.v_x)):
            for ni in range(-1, 2):
                for nj in range(-1, 2):
                    for j in range(len(self.x)):
                        dx = self.v_x[i] - self.x[j] + float(ni * self._Lx)
                        dy = self.v_y[i] - self.y[j] + float(nj * self._Ly)
                        dist = np.sqrt(dx ** 2 + dy ** 2)
                        if ((i + 1) % 3 != 0) and (1.1 > dist > 0.1):
                            dx /= dist
                            dy /= dist
                            lij = (dx * self.mx[j] + dy * self.my[j]) / dist ** 2
                            i_vertex.append(i)
                            j_vertex.append(j)
                            i_vertex_q.append(lij)
                        elif 1.3 > dist > 0.9:
                            dx /= dist
                            dy /= dist
                            lij = (dx * self.mx[j] + dy * self.my[j]) / dist ** 2
                            i_vertex.append(i)
                            j_vertex.append(j)
                            i_vertex_q.append(np.around(lij, 2))

        self.i_vertex = np.array(i_vertex, dtype=int)
        self.j_vertex = np.array(j_vertex, dtype=int)
        self.i_vertex_q = np.array(i_vertex_q)

    def _init_grid(self):
        xmin, xmax = min(self.x), max(self.x)
        ymin, ymax = min(self.y), max(self.y)
        x = np.linspace(xmin, xmax, self.nx)
        y = np.linspace(ymin, ymax, self.ny)
        self.XX, self.YY = np.meshgrid(x, y)

    def set_spin(self, spin):
        self.Spin = np.array(spin, dtype=int)
        self.calc_charge()

    def calc_charge(self):
        self.charge = np.zeros(shape=len(self.v_x))
        for i, j, q in zip(self.i_vertex, self.j_vertex, self.i_vertex_q):
            self.charge[i] += q * self.Spin[j]

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

    def plot_model(self, **kwargs):
        pos = np.column_stack([self.x, self.y])
        vec = np.column_stack([self.mx * self.Spin, self.my * self.Spin])
        return plot_vector(pos, vec, **kwargs)

    def plot_model_vertice(self, label=False, ax=None):
        self.calc_charge()
        return plot_vertice(self.v_x, self.v_y, self.charge, ax=ax, label=label)
