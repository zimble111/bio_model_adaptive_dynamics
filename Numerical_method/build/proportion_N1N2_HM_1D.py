import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("HM_1D.csv")

ratio = data['N 1'] / data['N 2']
new_data = data[np.isfinite(ratio)]
plt.title('Proportions of N1/N2 for HM 1D')
contour = plt.tricontourf(new_data['sd_in'], new_data['sd_ex'], ratio[np.isfinite(ratio)], cmap='Greens', levels=25)
plt.colorbar(contour)
plt.xlabel('sd_w_in')
plt.ylabel('sd_w_ex')
plt.plot([], [], ' ', label='Proportions of N1/N2')
plt.legend(loc=2)
plt.savefig("Proportions_of_N1_N2_for_HM_1D.png")
plt.show()
