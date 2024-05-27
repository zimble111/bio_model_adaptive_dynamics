import subprocess
from itertools import product
import numpy as np
import json
from copy import deepcopy
import cbor2
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly
import plotly.express as px
from os import path

COLORS = {
    'N1': 'b',
    'N2': 'r',
}

COOLORS_SCALE = (
    'Blues',
    'Reds',
)

'''
CMAPS = {
    'N1': cm.Blues,
    'N2': cm.Reds,
}
'''

LABELS_LATEX = {
    'N1': '$N_1$',
    'N2': '$N_2$',
}

def GetPlotParams(n):
    return {
        'label': LABELS_LATEX[n],
        'c': COLORS[n],
    }

def GetPlotParams3D(n):
    return {
        #'label': LABELS_LATEX[n],
        #'color': COLORS[n],
    }

FILE_PATH = 'CCTO_3D.json'

with open(FILE_PATH, 'r') as f:
    params = json.load(f)

base_shape = (params['x_len'], params['y_len'])
x_vec = np.ndarray(shape=base_shape[0])
y_vec = np.ndarray(shape=base_shape[1])

data_names = params.pop('datas')
for data_name in data_names:
    x_vec[data_name['i']] = data_name['x']
    y_vec[data_name['j']] = data_name['y']

points = params['points']
species = params['species']
title = params['title']
x_label = params['x_label']
y_label = params['y_label']

print(json.dumps(params, indent=4))
FOLDER = path.join(path.dirname(FILE_PATH), params['title'])

def get2DArray():
    return np.ndarray(shape=base_shape)


X, Y = np.meshgrid(x_vec, y_vec, indexing='ij')
Ns = np.ndarray(shape=(params['species'],) + base_shape)
Ds = np.ndarray(shape=(species, species,) + base_shape + (points,))
N1s = []
N2s = []
def FillDatas():
    for data_name in data_names:
        i = data_name['i']
        j = data_name['j']
        try:
            with open('{0}/{1}.N'.format(FOLDER, data_name['name']), 'r') as f:
                s = f.readline().split()
                N1 = 0 if s[0] == '-nan' else float(s[0])
                N2 = 0 if s[1] == '-nan' else float(s[1])
                N1 = max(N1, 0)
                N2 = max(N2, 0)
                Ns[0, i, j] = N1
                Ns[1, i, j] = N2
            with open('{0}/{1}.data'.format(FOLDER, data_name['name']), 'r') as f:
                Ds[0, 0, i, j, :] = np.array(f.readline().split(), dtype=float)
                Ds[0, 1, i, j, :] = np.array(f.readline().split(), dtype=float)
                Ds[1, 1, i, j, :] = np.array(f.readline().split(), dtype=float)
        except Exception as e:
            print(e)
            pass

print(Ds.shape)
print(Ns.shape)
FillDatas()

fig = go.Figure()
for i in range(species):
    fig.add_surface(z=Ns[i], x=X, y=Y, name='N{}'.format(i+1), colorscale=COOLORS_SCALE[i], showscale=False)
fig.update_layout(
    title='Ns for {title}'.format(title=title),
    autosize=True,
    width=600,
    height=600,
    margin=dict(l=65, r=50, b=65, t=90),
    scene={
        'xaxis_title': x_label,
        'yaxis_title': y_label,
        'zaxis_title': 'N',
    }
)
fig.show()

contour = np.ndarray(shape=base_shape + (3,), dtype=np.uint8)
eps = 0.001

contour = get2DArray()
contour[np.bitwise_and(Ns[0] < eps, Ns[1] < eps)] = 0     # Оба вымерли
contour[np.bitwise_and(Ns[0] >= eps, Ns[1] < eps)] = 1    # Выжил только первый вид
contour[np.bitwise_and(Ns[0] < eps, Ns[1] >= eps)] = 2    # Выжил только второй вид
contour[np.bitwise_and(Ns[0] >= eps, Ns[1] >= eps)] = 3   # Сосуществование
contour = np.rot90(np.flip(contour, (1)))

fig = go.Figure()
fig.add_heatmap(
    z=contour,
    x=x_vec,
    y=y_vec,
    zmin=0,
    zmax=3,
)
fig.update_layout(
    title='Ns for {title}'.format(title=title),
    autosize=True,
    width=800,
    height=800,
    margin=dict(l=65, r=50, b=65, t=90),
    scene={
        'xaxis_title': x_label,
        'yaxis_title': y_label,
        'zaxis_title': 'N',
    }
)
fig.show()