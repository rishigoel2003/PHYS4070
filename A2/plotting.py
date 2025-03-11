import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# Load the data
data = pd.read_csv('IsingModel.txt', sep=' ', names=['index', 'temperature', 'energy', 'magnetization'])

# Get unique temperatures
temperatures = np.sort(data['temperature'].unique())

# Initialize arrays to store results
avg_energy = []
sem_energy = []  # Standard error of the mean for energy
avg_mag = []
sem_mag = []     # Standard error of the mean for magnetization




# Process data for each temperature
for temp in temperatures:
    # Get data for this temperature
    temp_data = data[data['temperature'] == temp]
    
    # Calculate average and standard error
    avg_energy.append(np.mean(temp_data['energy']))
    sem_energy.append(stats.sem(temp_data['energy']))
    
    # For magnetization, use absolute values to handle phase symmetry
    mag_values = np.abs(temp_data['magnetization'])
    avg_mag.append(np.mean(mag_values))
    sem_mag.append(stats.sem(mag_values))






# Create figure with two subplots
plt.figure(figsize=(12, 5),num=1)





# Plot energy per spin vs temperature
plt.subplot(1, 2, 1)
plt.errorbar(temperatures, avg_energy, yerr=sem_energy, fmt='o-', capsize=3)
plt.xlabel('Temperature')
plt.ylabel('Energy per Spin')
plt.title('Energy vs Temperature')
plt.grid(True, alpha=0.3)




# Define Onsager's exact solution for magnetization
def onsager_magnetization(T, J=1, kb=1):
    # Theoretical critical temperature
    Tc = 2.269185 * J / kb
    
    # Magnetization is zero above Tc
    if T >= Tc:
        return 0
    
    # Onsager's formula for T < Tc
    return (1 - (np.sinh(2*J/(kb*T)))**(-4))**(1/8)

# Calculate theoretical magnetization for a range of temperatures
theory_temps = np.linspace(min(temperatures), max(temperatures), 100)
theory_mag = np.array([onsager_magnetization(T) for T in theory_temps])






# Plot magnetization per spin vs temperature with theory curve
plt.subplot(1, 2, 2)
plt.errorbar(temperatures, avg_mag, yerr=sem_mag, fmt='o', capsize=3, label='Simulation')
plt.plot(theory_temps, theory_mag, 'r-', label="Onsager's exact solution")
plt.xlabel('Temperature')
plt.ylabel('|Magnetization| per Spin')
plt.title('Magnetization vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Plots/2.1_A.png', dpi=300)
plt.show()






# Find the critical temperature based on maximum in specific heat
energy_variance = []
for temp in temperatures:
    temp_data = data[data['temperature'] == temp]
    energy_variance.append(np.var(temp_data['energy']))

specific_heat = np.array(energy_variance) *16**2 / (temperatures ** 2)
critical_temp_index = np.argmax(specific_heat) 
critical_temp = temperatures[critical_temp_index]





print(f"Estimated critical temperature: {critical_temp:.3f}")
print(f"Theoretical critical temperature: 2.269")

# Plot specific heat
plt.figure(figsize=(8, 5),num=1)
plt.plot(temperatures, specific_heat, 'o-')
plt.axvline(x=critical_temp, color='r', linestyle='--', 
            label=f'Estimated critical T ≈ {critical_temp:.3f}')
plt.axvline(x=2.269, color='g', linestyle='--', 
            label=f'Theoretical critical T = 2.269')
plt.xlabel('Temperature')
plt.ylabel('Specific Heat per spin')
plt.title('Specific Heat vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('Plots/2.1_B.png', dpi=300)
plt.show()





# Calculate and plot susceptibility
susceptibility = []
for temp in temperatures:
    temp_data = data[data['temperature'] == temp]
    mag_values = np.abs(temp_data['magnetization']) 
    susceptibility.append(np.var(mag_values) *16**2 / temp)



plt.figure(figsize=(8, 5),num=1)
plt.plot(temperatures, susceptibility, 'o-')
plt.axvline(x=critical_temp, color='r', linestyle='--', 
            label=f'Estimated critical T ≈ {critical_temp:.3f}')
plt.axvline(x=2.269, color='g', linestyle='--', 
            label=f'Theoretical critical T = 2.269')
plt.xlabel('Temperature')
plt.ylabel('Magnetic Susceptibility per spin')
plt.title('Susceptibility vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('Plots/2.1_C.png', dpi=300)
plt.show()
