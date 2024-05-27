import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("HM_3D.csv")

plt.title("coexistence of N1, HM 3D")
plt.scatter(data.loc[data['N 1'] == 0, 'sd_ex'], data.loc[data['N 1'] == 0, 'sd_in'], color='red', label='N 1 = 0')
plt.scatter(data.loc[data['N 1'] != 0, 'sd_ex'], data.loc[data['N 1'] != 0, 'sd_in'], color='blue', label='N 1 != 0')
plt.ylabel('sd_in')
plt.xlabel('sd_ex')
plt.legend(["Not coexistence", "Coexistence"], loc=2)
plt.savefig("coex_of_N1_HM_3D.png")
plt.show()