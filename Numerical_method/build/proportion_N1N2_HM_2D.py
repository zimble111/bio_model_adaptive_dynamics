import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("HM_2D.csv")

plt.title('Proportions of N1/N2 for HM 2D')
plt.tricontourf(data['sd_in'], data['sd_ex'], data['N 1']/data['N 2'], camp='Greens', label='N1/N2', levels=25)
plt.colorbar()
plt.xlabel('sd_w_in')
plt.ylabel('sd_w_ex')
plt.savefig("Proportions_of_N1_N2_for_HM_2D.png")
plt.legend()
plt.show()