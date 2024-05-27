import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data = pd.read_csv("HM_1D.csv")
new_data = data[np.isfinite(data['N 2'])]
plt.title('HM 1D, type 2')
contour = plt.tricontourf(new_data['sd_in'], new_data['sd_ex'], new_data['N 2'], cmap='Blues', levels=25)
plt.colorbar(contour)
plt.xlabel('sd_in')
plt.ylabel('sd_ex')
plt.plot([], [], ' ', label='Competition of N2')
plt.legend(loc=2)
plt.savefig("HM_1D_type_2.png")
plt.show()
