import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("HM_1D.csv")

plt.title("coexistence of N2, HM 1D")
plt.scatter(data.loc[data['N 2'] == 0, 'sd_ex'], data.loc[data['N 2'] == 0, 'sd_in'], color='red', label='N 2 = 0')
plt.scatter(data.loc[data['N 2'] != 0, 'sd_ex'], data.loc[data['N 2'] != 0, 'sd_in'], color='blue', label='N 2 != 0')
plt.ylabel('sd_in')
plt.xlabel('sd_ex')
plt.legend(["Not coexistence", "Coexistence"], loc=2)
plt.savefig("coex_of_N2_HM_1D.png")
plt.show()