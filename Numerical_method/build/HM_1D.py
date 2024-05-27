import csv
import numpy as np
from subprocess import run, call


with open('HM_1D.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    field = ["sd_in", "sd_ex", "N 1", "N 2"]
    writer.writerow(field)
    i = 0
    j = 0
    for sd_in in np.linspace(0.01, 0.22, num=10):
        for sd_ex in np.linspace(0.01, 0.22, num=10):
            name = 'HM_1D/HM_1D_{i}_{j}'.format(i=i, j=j)
            cmd = ['/Users/pavelgnilomedov/Desktop/w/build/a', '-t' + name, '-species 2', '-dim 1', '-al 0.4', '-points 1024', '-a 1', '-b', '0.4', '0.4', '-dvec', '0.2', '0.2', '-dmat 0.001', '0.001', '0.001', '0.001', '-sw', str(sd_in), str(sd_ex), str(sd_ex), str(sd_in), '-sm', '0.04', '0.06']
            call(cmd)
            j += 1
        i += 1
    i = 0
    j = 0
    for sd_in in np.linspace(0.01, 0.22, num=10):
        for sd_ex in np.linspace(0.01, 0.22, num=10):
            s = np.genfromtxt('/Users/pavelgnilomedov/Desktop/w/build/HM_1D/HM_1D_{i}_{j}.N'.format(i=i, j=j))
            N1 = 0 if s[0] == 'nan' else float(s[0])
            N2 = 0 if s[1] == 'nan' else float(s[1])
            N1 = max(N1, 0)
            N2 = max(N2, 0)
            writer.writerow([str(sd_in), str(sd_ex), str(10 * N1), str(10 * N2)])
            j += 1
        i += 1
