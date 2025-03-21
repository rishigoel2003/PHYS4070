import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# Load the data
data = pd.read_csv('IsingAppliedB.txt', sep=' ', names=['index', 'B', 'temperature', 'energy', 'magnetization'])

# Get unique temperatures and lengths
temperatures = np.sort(data['temperature'].unique())
mags = np.sort(data['B'].unique())

# Create figure with two subplots
plt.figure(figsize=(12, 5))

# Plot energy vs temperature for each lattice size
plt.subplot(1, 2, 1)
for mag in mags:
    # Initialize arrays for this specific length
    avg_energy = []
    sem_energy = []
    
    # Process data for each temperature with this length
    for temp in temperatures:
        # Get data for this temperature and length
        temp_data = data[(data['temperature'] == temp) & (data['B'] == mag)]
        avg_energy.append(np.mean(temp_data['energy']))
        sem_energy.append(stats.sem(temp_data['energy']))
        
    # Plot energy for this length
    plt.errorbar(temperatures, avg_energy, yerr=sem_energy, fmt='o-', capsize=3, 
                 label=f'B = {float(mag)}')

plt.xlabel('Temperature')
plt.ylabel('Energy per Spin')
plt.title('Energy vs Temperature')
plt.legend()
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

# Plot magnetization vs temperature for each lattice size
plt.subplot(1, 2, 2)
for mag in mags:
    # Initialize arrays for this specific length
    avg_mag = []
    sem_mag = []
    
    # Process data for each temperature with this length
    for temp in temperatures:
        # Get data for this temperature and length
        temp_data = data[(data['temperature'] == temp) & (data['B'] == mag)]
        
        if not temp_data.empty:
            # For magnetization, use absolute values to handle phase symmetry
            mag_values = np.abs(temp_data['magnetization'])
            avg_mag.append(np.mean(mag_values))
            sem_mag.append(stats.sem(mag_values))
        else:
            # Handle missing data points
            avg_mag.append(np.nan)
            sem_mag.append(np.nan)
    
    # Plot magnetization for this length
    plt.errorbar(temperatures, avg_mag, yerr=sem_mag, fmt='o-', capsize=3, 
                 label=f'B = {float(mag)}')

# Add the theoretical curve
plt.plot(theory_temps, theory_mag, 'k--', label="Onsager's exact solution")
plt.xlabel('Temperature')
plt.ylabel('|Magnetization| per Spin')
plt.title('Magnetization vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Plots/B2_3.png', dpi=300)
# plt.show()








# Plot specific heat for each lattice size
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)

for B in mags:
    specific_heat = []
    
    for temp in temperatures:
        #take data
        temp_data = data[(data['temperature'] == temp) & (data['B'] == B)]
        energy_variance = np.var(temp_data['energy'])
        specific_heat.append(energy_variance *16**2 / (temp ** 2))


    plt.plot(temperatures, specific_heat, 'o-', label=f'B = {float(B)}')

plt.axvline(x=2.269, color='k', linestyle='--', 
            label=f'Theoretical critical T = 2.269')
plt.xlabel('Temperature')
plt.ylabel('Specific Heat')
plt.title('Specific Heat vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)




# Plot susceptibility for each lattice size
plt.subplot(1, 2, 2)

for B in mags:
    susceptibility = []
    
    for temp in temperatures:
        temp_data = data[(data['temperature'] == temp) & (data['B'] == B)]
        
        if not temp_data.empty:
            mag_values = np.abs(temp_data['magnetization'])
            mag_variance = np.var(mag_values)
            susceptibility.append(mag_variance *16**2  / temp)
        else:
            susceptibility.append(np.nan)
    
    plt.plot(temperatures, susceptibility, 'o-', label=f'B = {float(B)}')

#vertical line
plt.axvline(x=2.269, color='k', linestyle='--', 
            label=f'Theoretical critical T = 2.269')
plt.xlabel('Temperature')
plt.ylabel('Magnetic Susceptibility')
plt.title('Susceptibility vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Plots/B2_4.png', dpi=300)
# plt.show()




