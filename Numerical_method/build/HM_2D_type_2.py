import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("HM_2D.csv")

plt.title('HM 2D, type 2')
contour = plt.tricontourf(data['sd_in'], data['sd_ex'], data['N 2'], cmap='Blues', levels=25)
plt.colorbar(contour)
plt.xlabel('sd_in')
plt.ylabel('sd_ex')
plt.plot([], [], ' ', label='Competition of N2')
plt.legend(loc=2)
plt.savefig("HM_2D_type_2.png")
plt.show()
