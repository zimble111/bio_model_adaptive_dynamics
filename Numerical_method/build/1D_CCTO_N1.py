import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


data = pd.read_csv('CCTO_1D.csv')

clean_data = data.dropna(subset=['N 1'])
clean_data = clean_data[np.isfinite(clean_data['N 1'])]

plt.title("one-dimensional CCTO case (Number of the first type)")
contour = plt.tricontourf(clean_data['sd_m_2'], clean_data['d12'], clean_data['N 1'], cmap='Greens')
plt.colorbar(contour)

plt.ylabel('d12')
plt.xlabel('sd_m_2')

plt.plot([], [], ' ', label='Competition of N1')
plt.legend(loc=2)

plt.savefig("1d_CCTO_Case_N1.png")

plt.show()
