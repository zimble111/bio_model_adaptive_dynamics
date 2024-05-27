import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv('coex_N1N2.csv')
data['d12'] = data['d12'].replace(',', '.', regex=True).astype(float)

plt.title("graph of coexistence with population proportion N1/N2")
plt.scatter(data['d12'], data['sd_m_2'], c=data['N 1']/data['N 2'], cmap='Blues', label='N1/N2')
plt.colorbar()
plt.xlabel('d12')
plt.ylabel('sd_m_2')
plt.legend(loc=2)
plt.savefig("GraphOfCoexistenceWithPopulationProportion_N1_N2.png")
plt.show()