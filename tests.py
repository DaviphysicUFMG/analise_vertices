import numpy as np
import pandas as pd
from main import VerticeAnalysis
import matplotlib.pyplot as plt


config0 = pd.read_csv('config0.csv', sep="\s+")
x = config0.x.to_numpy()
y = config0.y.to_numpy()
mx = config0.mx.to_numpy()
my = config0.my.to_numpy()

model = VerticeAnalysis(x, y, mx, my)

carga_k = 0
carga_T = 0
for i in range(model.Ns):
    if (i + 1) % 3 != 0:
        carga_k += model.charge[i]
    else:
        carga_T += model.charge[i]

print("Carga Kagome:", np.around(carga_k, 3))
print("Carga Triangular:", np.around(carga_T, 3))

plt.figure(figsize=(12, 10))
ax = plt.gca()

scatter = model.plot_model_vertice(ax=ax, label=True)
quiv = model.plot_model(ax=ax, **{"cmap": 'peem30'})

plt.show()
