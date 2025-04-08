import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the data
data = pd.read_csv('Hydrogen_Like_Energy.txt', sep=' ', names=['l', 'n', 'energy'])

# Extract relevant columns
ls = data['l']
ns = data['n']
energies = data['energy']

# Select the first 10 values of n
first_10_data = data[data['n'] <= 10]

# Theoretical energy for n values
n_values = np.arange(1, 11)
theoretical_energies = 9 / (2 * n_values**2)

# Plot the data
plt.figure(figsize=(10, 5))

# Loop through unique values of l and plot the experimental energies for each l
unique_ls = first_10_data['l'].unique()
for l in unique_ls:
    l_data = first_10_data[first_10_data['l'] == l]
    plt.plot(l_data['n'], l_data['energy'], label=f'Experimental Energy (l={l})', marker='o')

# Plot the theoretical energy
plt.plot(n_values, theoretical_energies, color='red', label='Theoretical Energy (9/2n^2)', linestyle='--')

# Labels and title
plt.xlabel('Quantum number n')
plt.ylabel('Energy')
plt.title('Experimental vs. Theoretical Energy for n = 1 to 10 (Different l Values)')
plt.legend()

# Show plot
plt.grid(True)
plt.savefig('Plots/B1_1', dpi=300)
# plt.show()




























# Load the data
data = pd.read_csv('Hydrogen_Like_Radial.txt', sep=' ', names=['l', 'n', 'r_expectation'])

# Select the first 10 values of n
first_10_data = data[data['n'] <= 10]

# Atomic number for the hydrogen-like atom
Z = 3

# Plot the data
plt.figure(figsize=(10, 5))

# Loop through unique values of l and plot the calculated expectation values for each l
unique_ls = first_10_data['l'].unique()
for l in unique_ls:
    l_data = first_10_data[first_10_data['l'] == l]
    
    # Calculate theoretical expectation values for this l
    theoretical_r = (3 * l_data['n']**2 - l * (l + 1)) / (2 * Z)
    
    plt.plot(l_data['n'], l_data['r_expectation'], label=f'Calculated <r> (l={l})', marker='o')
    plt.plot(l_data['n'], theoretical_r, label=f'Theoretical <r> (l={l})', linestyle='--')

# Labels and title
plt.xlabel('Quantum number n')
plt.ylabel('Expectation value <r>')
plt.title('Calculated vs. Theoretical Radial Expectation Values for n = 1 to 10')
plt.legend()

# Show plot
plt.grid(True)
plt.savefig('Plots/B1_2', dpi=300)
# plt.show()



















import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the radial probability density data
data = pd.read_csv('Hydrogen_Like_Probability_Density.txt', sep=' ', 
                   names=['l', 'n', 'r', 'prob_density'])

# Create a new figure with a reasonable size
plt.figure(figsize=(10, 6))

# Plot for l=0, n=1 (ground state)
state_data = data[(data['l'] == 0) & (data['n'] == 1)]
plt.plot(state_data['r'], state_data['prob_density'], 
         label='l=0, n=1 (1s)', color='blue', linewidth=2)

# Plot for l=0, n=2 (first excited s state)
state_data = data[(data['l'] == 0) & (data['n'] == 2)]
plt.plot(state_data['r'], state_data['prob_density'], 
         label='l=0, n=2 (2s)', color='red', linewidth=2)

# Plot for l=1, n=2 (first p state)
state_data = data[(data['l'] == 1) & (data['n'] == 2)]
plt.plot(state_data['r'], state_data['prob_density'], 
         label='l=1, n=2 (2p)', color='green', linewidth=2)

# Add labels and title
plt.xlabel('Radial Distance (r)', fontsize=12)
plt.ylabel('Radial Probability Density', fontsize=12)
plt.title('Hydrogen-like Atom Radial Probability Densities (Z=3)', fontsize=14)

# Add legend
plt.legend(fontsize=10)

# Set reasonable x-axis limit (adjust based on your data)
plt.xlim(0, 15)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Save the figure
plt.savefig('Plots/B1_3', dpi=300)
# plt.show()
