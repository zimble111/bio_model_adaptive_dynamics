#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 1, 1024)
D = np.genfromtxt('HM_3D/HM_3D_4_42.data')
plt.plot(x, D[0], label=r'$D_{1, 1}$')
plt.plot(x, D[1], label=r'$D_{1, 2}$')
plt.plot(x, D[2], label=r'$D_{2, 2}$')
plt.legend()
plt.show()
