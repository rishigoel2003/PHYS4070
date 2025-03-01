import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

# Load the data
data = pd.read_csv('IsingPower.txt', sep=' ', names=['index', 'temperature', 'energy', 'magnetization'])

# Get unique temperatures and lengths
temperatures = np.sort(data['temperature'].unique())

avg_mag = []
sem_mag = []

# Process data for each temperature with this length
for temp in temperatures:
    # Get data for this temperature and length
    temp_data = data[(data['temperature'] == temp)]
    
    if not temp_data.empty:
        # For magnetization, use absolute values to handle phase symmetry
        mag_values = np.abs(temp_data['magnetization'])
        avg_mag.append(np.mean(mag_values))
        sem_mag.append(stats.sem(mag_values))

# Convert to numpy arrays
temperatures = np.array(temperatures)
avg_mag = np.array(avg_mag)
sem_mag = np.array(sem_mag)
temp_diff = np.abs(temperatures - 2.269)  # Distance from critical temperature

# Calculate susceptibility
susceptibility = []
for temp in temperatures:
    temp_data = data[data['temperature'] == temp]
    mag_values = np.abs(temp_data['magnetization']) 
    susceptibility.append(np.var(mag_values) *16**2 / temp)

susceptibility = np.array(susceptibility)

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Filter for magnetization plot (T < Tc)
valid_mag_indices = (avg_mag > 0) & (temp_diff > 0) & (temperatures < 2.269)
x_mag = np.log10(temp_diff[valid_mag_indices])
y_mag = np.log10(avg_mag[valid_mag_indices])

# Perform linear regression for magnetization
slope_mag, intercept_mag, r_value_mag, p_value_mag, std_err_mag = stats.linregress(x_mag, y_mag)

# Theoretical values for 2D Ising model
beta_theory = 1/8  # = 0.125
gamma_theory = 7/4  # = 1.75

# Plot magnetization data and regression line
ax1.scatter(temp_diff[valid_mag_indices], avg_mag[valid_mag_indices], marker='o', label='Simulation Data')
x_line = np.logspace(np.log10(min(temp_diff[valid_mag_indices])), np.log10(max(temp_diff[valid_mag_indices])), 100)
y_line = 10**intercept_mag * x_line**slope_mag
ax1.loglog(x_line, y_line, 'r-', label=f'Fit: β = {slope_mag:.3f}±{std_err_mag:.3f}')

# Add theoretical line for magnetization
# Use the same amplitude (10**intercept_mag) but with theoretical exponent
y_theory = 10**intercept_mag * x_line**beta_theory
ax1.loglog(x_line, y_theory, 'g--', label=f'Theory: β = {beta_theory}')

ax1.set_xlabel('|T - T_c|')
ax1.set_ylabel('|Magnetization| per Spin')
ax1.set_title('Magnetization vs |T - T_c| (Log-Log)')
ax1.grid(True, alpha=0.3, which="both")
ax1.legend()

# Filter for susceptibility plot
valid_sus_indices = (susceptibility > 0) & (temp_diff > 0)
x_sus = np.log10(temp_diff[valid_sus_indices])
y_sus = np.log10(susceptibility[valid_sus_indices])

# Perform linear regression for susceptibility
slope_sus, intercept_sus, r_value_sus, p_value_sus, std_err_sus = stats.linregress(x_sus, y_sus)

# Plot susceptibility data and regression line
ax2.scatter(temp_diff[valid_sus_indices], susceptibility[valid_sus_indices], marker='o', label='Simulation Data')
x_line = np.logspace(np.log10(min(temp_diff[valid_sus_indices])), np.log10(max(temp_diff[valid_sus_indices])), 100)
y_line = 10**intercept_sus * x_line**slope_sus
ax2.loglog(x_line, y_line, 'r-', label=f'Fit: γ = {-slope_sus:.3f}±{std_err_sus:.3f}')

# Add theoretical line for susceptibility
# Use the same amplitude (10**intercept_sus) but with theoretical exponent
y_theory = 10**intercept_sus * x_line**(-gamma_theory)
ax2.loglog(x_line, y_theory, 'g--', label=f'Theory: γ = {gamma_theory}')

ax2.set_xlabel('|T - T_c|')
ax2.set_ylabel('Magnetic Susceptibility')
ax2.set_title('Susceptibility vs |T - T_c| (Log-Log)')
ax2.grid(True, alpha=0.3, which="both")
ax2.legend()

plt.tight_layout()
# plt.savefig('ising_critical_exponents_with_theory.png', dpi=300)
plt.show()

# Print the results and comparison with theory
print(f"Magnetization critical exponent (β):")
print(f"  Measured: {slope_mag:.4f} ± {std_err_mag:.4f}")
print(f"  Theoretical: {beta_theory:.4f}")
print(f"  Difference: {abs(slope_mag - beta_theory):.4f}")
print(f"  R-squared: {r_value_mag**2:.4f}")

print(f"\nSusceptibility critical exponent (γ):")
print(f"  Measured: {-slope_sus:.4f} ± {std_err_sus:.4f}")
print(f"  Theoretical: {gamma_theory:.4f}")
print(f"  Difference: {abs(-slope_sus - gamma_theory):.4f}")
print(f"  R-squared: {r_value_sus**2:.4f}")