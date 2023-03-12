import numpy as np
import pandas as pd
from main import VerticeAnalysis
from plot import plot_vertice
import matplotlib.pyplot as plt


config0 = pd.read_csv('config0.csv')
x = config0.x.to_numpy()
y = config0.y.to_numpy()
mx = config0.mx.to_numpy()
my = config0.my.to_numpy()

model = VerticeAnalysis(x, y, mx, my, nx=100, ny=100)

XX = model.XX
YY = model.YY
Bx0 = model.Bx.copy()
By0 = model.By.copy()
charge0 = model.charge.copy()

model.Spin[133] *= -1
model.Spin[162] *= -1
model.Spin[163] *= -1
model.Spin[195] *= -1
model.Spin[196] *= -1

model.calc_field()
Bx = model.Bx
By = model.By

Bx1 = Bx - Bx0
By1 = By - By0

kwargs = {
    'cmap': 'flatspin',
    'alpha': 1.,
}
model.get_charge()
charge1 = model.charge - charge0

fig, axes = plt.subplots(1, 2, figsize=(8, 3))

quiv = model.plot_model(ax=axes[0], **kwargs)
plot_vertice(model.v_x, model.v_y, charge1, ax=axes[1])
axes[1].streamplot(XX, YY, Bx1, By1, density=1, linewidth=0.8)
axes[1].axis('equal')
axes[1].axis('off')
plt.show()


# carga_k = 0
# carga_T = 0
# for i in range(model.Ns):
#     if (i + 1) % 3 != 0:
#         carga_k += model.charge[i]
#     else:
#         carga_T += model.charge[i]
#
# print("Carga Kagome:", np.around(carga_k, 3))
# print("Carga Triangular:", np.around(carga_T, 3))
#
# plt.figure(figsize=(12, 10))
# ax = plt.gca()
#
# scatter = model.plot_model_vertice(ax=ax, label=True)
# quiv = model.plot_model(ax=ax, **{"cmap": 'peem30'})
#
# plt.show()
