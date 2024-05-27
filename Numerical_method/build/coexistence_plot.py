import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv('coex_N1N2.csv')
data['d12'] = data['d12'].replace(',', '.', regex=True).astype(float)

plt.title("graph of coexistence")
plt.scatter(data.loc[data['N 1'] == 0, 'd12'], data.loc[data['N 1'] == 0, 'sd_m_2'], color='black', label='N 1 = 0')
plt.scatter(data.loc[(data['N 1'] != 0) & (data['N 2'] != 0), 'd12'], data.loc[(data['N 1'] != 0) & (data['N 2'] != 0), 'sd_m_2'], color='yellow', label='N1 != 0 && N2 !=0')
plt.scatter(data.loc[data['N 2'] == 0, 'd12'], data.loc[data['N 2'] == 0, 'sd_m_2'], color='red', label='N 2 = 0')
plt.ylabel('sd_m_2')
plt.xlabel('d12')
plt.legend(["Victim","Coexistence", "Predator"], loc=1)
plt.show()