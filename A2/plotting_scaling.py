import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# Load the data
data = pd.read_csv('IsingScaling.txt', sep=' ', names=['index', 'length', 'temperature', 'energy', 'magnetization'])

# Get unique temperatures and lengths
temperatures = np.sort(data['temperature'].unique())
lengths = np.sort(data['length'].unique())

# Create figure with two subplots
plt.figure(figsize=(12, 5))

# Plot energy vs temperature for each lattice size
plt.subplot(1, 2, 1)
for length in lengths:
    # Initialize arrays for this specific length
    avg_energy = []
    sem_energy = []
    
    # Process data for each temperature with this length
    for temp in temperatures:
        # Get data for this temperature and length
        temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
        if not temp_data.empty:
            # Calculate average and standard error
            avg_energy.append(np.mean(temp_data['energy']))
            sem_energy.append(stats.sem(temp_data['energy']))
        else:
            # Handle missing data points
            avg_energy.append(np.nan)
            sem_energy.append(np.nan)
    
    # Plot energy for this length
    plt.errorbar(temperatures, avg_energy, yerr=sem_energy, fmt='o-', capsize=3, 
                 label=f'L = {int(length)}')

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
for length in lengths:
    # Initialize arrays for this specific length
    avg_mag = []
    sem_mag = []
    
    # Process data for each temperature with this length
    for temp in temperatures:
        # Get data for this temperature and length
        temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
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
                 label=f'L = {int(length)}')

# Add the theoretical curve
plt.plot(theory_temps, theory_mag, 'k--', label="Onsager's exact solution")
plt.xlabel('Temperature')
plt.ylabel('|Magnetization| per Spin')
plt.title('Magnetization vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
# plt.savefig('ising_model_results.png', dpi=300)
plt.show()








# # Plot specific heat for each lattice size
# plt.figure(figsize=(12, 5))
# plt.subplot(1, 2, 1)

# for length in lengths:
#     specific_heat = []
    
#     for temp in temperatures:
#         temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
#         if not temp_data.empty:
#             energy_variance = np.var(temp_data['energy'])
#             specific_heat.append(energy_variance *length**2 / (temp ** 2))
#         else:
#             specific_heat.append(np.nan)
    
#     plt.plot(temperatures, specific_heat, 'o-', label=f'L = {int(length)}')

# plt.axvline(x=2.269, color='k', linestyle='--', 
#             label=f'Theoretical critical T = 2.269')
# plt.xlabel('Temperature')
# plt.ylabel('Specific Heat')
# plt.title('Specific Heat vs Temperature')
# plt.legend()
# plt.grid(True, alpha=0.3)

# # Plot susceptibility for each lattice size
# plt.subplot(1, 2, 2)

# for length in lengths:
#     susceptibility = []
    
#     for temp in temperatures:
#         temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
#         if not temp_data.empty:
#             mag_values = np.abs(temp_data['magnetization'])
#             mag_variance = np.var(mag_values)
#             susceptibility.append(mag_variance *length**2  / temp)
#         else:
#             susceptibility.append(np.nan)
    
#     plt.plot(temperatures, susceptibility, 'o-', label=f'L = {int(length)}')

# plt.axvline(x=2.269, color='k', linestyle='--', 
#             label=f'Theoretical critical T = 2.269')
# plt.xlabel('Temperature')
# plt.ylabel('Magnetic Susceptibility')
# plt.title('Susceptibility vs Temperature')
# plt.legend()
# plt.grid(True, alpha=0.3)

# plt.tight_layout()
# # plt.savefig('ising_model_thermodynamics.png', dpi=300)
# plt.show()

















# Plot specific heat and susceptibility with finite-size scaling
plt.figure(figsize=(12, 10))

# Critical temperature
Tc = 2.269

# Critical exponents for 2D Ising model
nu = 1
gamma = 7/4

# First row: Original plots
# Specific Heat
plt.subplot(2, 2, 1)
for length in lengths:
    specific_heat = []
    
    for temp in temperatures:
        temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
        if not temp_data.empty:
            energy_variance = np.var(temp_data['energy'])
            specific_heat.append(energy_variance * length**2 / (temp ** 2))
        else:
            specific_heat.append(np.nan)
    
    plt.plot(temperatures, specific_heat, 'o-', label=f'L = {int(length)}')

plt.axvline(x=Tc, color='k', linestyle='--', label=f'Critical T = {Tc}')
plt.xlabel('Temperature')
plt.ylabel('Specific Heat')
plt.title('Original: Specific Heat vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)

# Susceptibility
plt.subplot(2, 2, 2)
for length in lengths:
    susceptibility = []
    
    for temp in temperatures:
        temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
        if not temp_data.empty:
            mag_values = np.abs(temp_data['magnetization'])
            mag_variance = np.var(mag_values)
            susceptibility.append(mag_variance * length**2 / temp)
        else:
            susceptibility.append(np.nan)
    
    plt.plot(temperatures, susceptibility, 'o-', label=f'L = {int(length)}')

plt.axvline(x=Tc, color='k', linestyle='--', label=f'Critical T = {Tc}')
plt.xlabel('Temperature')
plt.ylabel('Magnetic Susceptibility')
plt.title('Original: Susceptibility vs Temperature')
plt.legend()
plt.grid(True, alpha=0.3)

# Second row: Rescaled plots
# Rescaled Specific Heat
plt.subplot(2, 2, 3)
for length in lengths:
    specific_heat = []
    reduced_temp = []
    
    for temp in temperatures:
        temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
        if not temp_data.empty:
            energy_variance = np.var(temp_data['energy'])
            c = energy_variance * length**2 / (temp ** 2)
            
            # Compute reduced temperature t = (T-Tc)/Tc
            t = (temp - Tc) / Tc
            
            specific_heat.append(c / np.log(length))  # c(t) = log(L)f(L^(1/ν)t)
            reduced_temp.append(length**(1/nu) * t)   # x-axis: L^(1/ν)t
        else:
            specific_heat.append(np.nan)
            reduced_temp.append(np.nan)
    
    plt.plot(reduced_temp, specific_heat, 'o-', label=f'L = {int(length)}')

plt.axvline(x=0, color='k', linestyle='--', label='Critical point')
plt.xlabel(r'Scaled Temperature: $L^{1/\nu}(T-T_c)/T_c$')
plt.ylabel(r'Scaled Specific Heat: $c(t)/\log(L)$')
plt.title('Finite-Size Scaling: Specific Heat')
plt.legend()
plt.grid(True, alpha=0.3)

# Rescaled Susceptibility
plt.subplot(2, 2, 4)
for length in lengths:
    susceptibility = []
    reduced_temp = []
    
    for temp in temperatures:
        temp_data = data[(data['temperature'] == temp) & (data['length'] == length)]
        
        if not temp_data.empty:
            mag_values = np.abs(temp_data['magnetization'])
            mag_variance = np.var(mag_values)
            chi = mag_variance * length**2 / temp
            
            # Compute reduced temperature t = (T-Tc)/Tc
            t = (temp - Tc) / Tc
            
            susceptibility.append(chi / length**gamma)  # χ(t) = L^γ/ν g(L^(1/ν)t)
            reduced_temp.append(length**(1/nu) * t)     # x-axis: L^(1/ν)t
        else:
            susceptibility.append(np.nan)
            reduced_temp.append(np.nan)
    
    plt.plot(reduced_temp, susceptibility, 'o-', label=f'L = {int(length)}')

plt.axvline(x=0, color='k', linestyle='--', label='Critical point')
plt.xlabel(r'Scaled Temperature: $L^{1/\nu}(T-T_c)/T_c$')
plt.ylabel(r'Scaled Susceptibility: $\chi(t)/L^{\gamma/\nu}$')
plt.title('Finite-Size Scaling: Susceptibility')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('ising_model_finite_size_scaling.png', dpi=300)
plt.show()