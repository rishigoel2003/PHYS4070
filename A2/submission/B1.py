import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the data
data = pd.read_csv('IsingB1.txt', sep=' ', names=['index', 'temperature', 'energy', 'magnetization'])

# Extract relevant columns
sweep_number = data['index']
energy = data['energy']
magnetization = data['magnetization']


# Compute cumulative mean
cumulative_mean_energy = np.cumsum(energy) / np.arange(1, len(energy) + 1)
cumulative_mean_magnetization = np.cumsum(magnetization) / np.arange(1, len(magnetization) + 1)

# Plot energy and magnetization
plt.figure(figsize=(10, 5))

# Energy plot
plt.subplot(1, 2, 1)
plt.plot(sweep_number, energy, label='Energy', alpha=0.5)
plt.plot(sweep_number, cumulative_mean_energy, label='Cumulative Mean Energy', linestyle='dashed', color='red')
plt.xscale('log')
plt.xlabel('Sweep Number')
plt.ylabel('Energy')
plt.legend()



# Magnetization plot
plt.subplot(1, 2, 2)
plt.plot(sweep_number, magnetization, label='Magnetization', alpha=0.5)
plt.plot(sweep_number, cumulative_mean_magnetization, label='Cumulative Mean Magnetization', linestyle='dashed', color='red')
plt.xscale('log')
plt.xlabel('Sweep Number')
plt.ylabel('Magnetization')
plt.legend()

plt.tight_layout()
plt.savefig('Plots/B1.png', dpi=300)
# plt.show()
