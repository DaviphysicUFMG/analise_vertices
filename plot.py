import numpy as np
from scipy.spatial import cKDTree
from numpy.linalg import norm
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.cm

fs_colours = [[.9, .1, .1],
              [.9, .9, .1],
              [.1, .7, .2],
              [.1, .7, .9],
              [.1, .1, .9],
              [.9, .1, .9],
              [.9, .1, .1]]
matplotlib.cm.register_cmap('flatspin', LinearSegmentedColormap.from_list('flatspin', fs_colours, N=1000))
normalize_rad = Normalize(vmin=0, vmax=2 * np.pi)

# Colormap from Li. et al., (2019)
# https://pubs.acs.org/doi/10.1021/acsnano.8b08884

clist = list('ckkbbkkggkkrrkkc')
clist = ['#00da00' if c == 'g' else c for c in clist]
clist = ['#0800da' if c == 'b' else c for c in clist]
clist = ['#ed0912' if c == 'r' else c for c in clist]
clist = ['#00ccff' if c == 'c' else c for c in clist]

cmap = ListedColormap(clist)
matplotlib.cm.unregister_cmap('li2019')
matplotlib.cm.register_cmap('li2019', cmap)


def vector_colors(U, V):
    C = np.arctan2(V, U)
    C[C < 0] = 2 * np.pi + C[C < 0]
    return C


def plot_vector(XY, UV, C=None, normalize=False, ax=None, **kwargs):
    XY = np.atleast_2d(XY)
    UV = np.atleast_2d(UV)
    XY = XY.reshape((-1, 2))
    X = XY[..., 0]
    Y = XY[..., 1]

    if normalize:
        nmax = norm(UV.reshape((-1, 2)), axis=-1).max(initial=None)
        if nmax != 0:
            UV = UV / nmax

    U = UV[..., 0]
    V = UV[..., 1]
    if C is None:
        C = vector_colors(U, V)
        kwargs.setdefault('norm', normalize_rad)

    if ax is None:
        ax = plt.gca()

    tree = cKDTree(XY)
    min_dist = np.mean(tree.query(XY, k=[2])[0])
    scale = 1.15 / min_dist if not np.isinf(min_dist) and min_dist != 0 else 1.0
    kwargs.setdefault('scale', scale)
    kwargs.setdefault('pivot', 'middle')
    kwargs.setdefault('width', 0.2 / scale)
    kwargs.setdefault('headwidth', 3)
    kwargs.setdefault('headlength', 2)
    kwargs.setdefault('headaxislength', 2)

    if type(kwargs["cmap"]) is str and kwargs["cmap"].startswith("peem"):
        peem_angle = np.deg2rad(float(kwargs["cmap"].strip("peem")))
        intensity = np.cos(np.linspace(0, 2 * np.pi, 360) - peem_angle)
        intensity = intensity * 0.5 + 0.5
        intensity = np.tile(intensity, (3, 1)).T
        kwargs["cmap"] = ListedColormap(intensity)
        ax.set_facecolor([0.5] * 3)

    quiv = ax.quiver(X, Y, U, V, C, units='xy',
                     angles='xy', scale_units='xy', **kwargs)
    ax.set_aspect('equal')

    xmin, xmax = np.min(X), np.max(X)
    ymin, ymax = np.min(Y), np.max(Y)

    ax.set_xlim(xmin - 1 / scale, xmax + 1 / scale)
    ax.set_ylim(ymin - 1 / scale, ymax + 1 / scale)

    return quiv


def plot_vertice(v_x, v_y, charge, ax=None, label=False, **kwargs):
    kwargs.setdefault('norm', Normalize(vmin=-1, vmax=1))
    kwargs.setdefault('cmap', 'bwr')

    if ax is None:
        ax = plt.gca()

    scatter = ax.scatter(v_x, v_y, c=charge, **kwargs)
    if label:
        for i in range(len(charge)):
            ax.text(v_x[i]+0.3, v_y[i], s=np.around(charge[i], 1))

    ax.set_aspect('equal')
    return scatter
