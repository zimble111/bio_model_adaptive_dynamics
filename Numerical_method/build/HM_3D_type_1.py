import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data = pd.read_csv("HM_3D.csv")
new_data = data[np.isfinite(data['N 1'])]
plt.title('HM 3D, type 1')
contour = plt.tricontourf(new_data['sd_in'], new_data['sd_ex'], new_data['N 1'], cmap='Blues', levels=25)
plt.colorbar(contour)
plt.xlabel('sd_in')
plt.ylabel('sd_ex')
plt.plot([], [], ' ', label='Competition of N1')
plt.legend(loc=2)
plt.savefig("HM_3D_type_1.png")
plt.show()
