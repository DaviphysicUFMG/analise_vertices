import numpy as np
import pandas as pd
from main import VerticeAnalysis
from plot import plot_vertice
import matplotlib.pyplot as plt

config0 = pd.read_csv('dados/config0.csv')
x = config0.x.to_numpy()
y = config0.y.to_numpy()
mx = config0.mx.to_numpy()
my = config0.my.to_numpy()

model = VerticeAnalysis(x, y, mx, my)

kwargs = {
    'cmap': 'peem45',
    'alpha': 1.,
}

spin = np.ones(len(x))
spin[0] = -1
spin[2] = -1
spin[4] = -1
model.set_spin(spin)

plt.figure(figsize=(12, 10))

model.plot_model(**kwargs)
plot_vertice(model.v_x, model.v_y, model.charge/3, label=True)
plt.show()
# Bx_init = model.Bx.copy()   # Campo magnético inicial
# By_init = model.By.copy()
# charge_init = model.charge.copy()   # Carga magnética inicial
#
# spins = pd.read_parquet('dados/Spin.parquet')
# spins = spins.T
#
# model.Spin = spins[0].to_list()
# model.calc_field()
# model.get_charge()
#
# Bx = model.Bx - Bx_init
# By = model.By - By_init
#
# charge = model.charge - charge_init
#
# kwargs = {
#     'cmap': 'flatspin',
#     'alpha': 1.,
# }
#
# fig, axes = plt.subplots(1, 2, figsize=(8, 3))
#
# quiv = model.plot_model(ax=axes[0], **kwargs)
# plot_vertice(model.v_x, model.v_y, charge, ax=axes[1])
# axes[1].streamplot(model.XX, model.YY, Bx, By, density=1, linewidth=0.8)
# axes[1].axis('equal')
# axes[1].axis('off')
# plt.show()
