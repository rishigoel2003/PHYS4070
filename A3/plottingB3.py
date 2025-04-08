import pandas as pd
import numpy as np
import matplotlib.pyplot as plt





# Load data
data = pd.read_csv('Hartree_Energy.txt', sep=' ', names=['iters','l', 'n', 'energy'])

# Filter data for n=1, l=0
n1_l0_data = data[(data['n'] == 1) & (data['l'] == 0)]

# Plot the data
plt.figure(figsize=(10, 6))

# Plot energy vs iterations for n=1, l=0
plt.plot(n1_l0_data['iters'], n1_l0_data['energy'], marker='o', linestyle='-', color='red')

# Labels and title
plt.xlabel('Iterations')
plt.ylabel('Energy (au)')
plt.title('Energy Convergence for 1s Orbital (n=1, l=0)')
# plt.legend()

# Improve readability
plt.grid(True)
plt.tight_layout()

# Save figure
plt.savefig('Plots/B3_1', dpi=300)
# plt.show()







# Load Hartree energy data
data = pd.read_csv('Hartree_Energy.txt', sep=' ', names=['iters', 'l', 'n', 'energy'])

# Get max iteration
max_iteration = data['iters'].max()

# Filter data for the final iteration and specific orbitals
final_iteration_data = data[(data['iters'] == max_iteration) & ((data['n'] == 2) & ((data['l'] == 0) | (data['l'] == 1)))]

# Theoretical energy values for 2s and 2p in atomic units
theoretical_energies = {'2s': -0.19814, '2p': -0.13023}

# Create figure
plt.figure(figsize=(10, 6))

# Scatter plot of final energies
colors = {0: 'green', 1: 'blue'}
for l_val in [0, 1]:
    subset = final_iteration_data[final_iteration_data['l'] == l_val]
    plt.scatter(
        subset['n'], 
        subset['energy'], 
        color=colors[l_val], 
        s=100, 
        label=f'Computed Energy (l={l_val})'
    )

# Plot theoretical values for 2s and 2p
plt.axhline(y=theoretical_energies['2s'], color='green', linestyle='--', label='Experimental 2s Energy')
plt.axhline(y=theoretical_energies['2p'], color='blue', linestyle='--', label='Experimental 2p Energy')

# Labels and title
plt.xlabel('Orbital (n)')
plt.ylabel('Energy (au)')
plt.title(f'Orbital Energies at Final Iteration (Iteration {max_iteration})')

# Add legend
plt.legend()

# Improve readability
plt.grid(True)
plt.tight_layout()

# Save figure
plt.savefig('Plots/B3_2', dpi=300)
# plt.show()










































# Load the radial probability density data
data = pd.read_csv('Hartree_Probability_Density.txt', sep=' ', 
                   names=['iters','l', 'n', 'r', 'prob_density'])

number_iters = data["iters"].unique().max()


data = data[data['iters'] == number_iters]

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
plt.title('Hartree Radial Probability Densities', fontsize=14)

# Add legend
plt.legend(fontsize=10)

plt.xlim(0, 15)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Save the figure
plt.savefig('Plots/B3_3', dpi=300)
# plt.show()









# Load the radial probability density data
data = pd.read_csv('Hartree_Probability_Density.txt', sep=' ', 
                   names=['iters','l', 'n', 'r', 'prob_density'])

# Filter data for n=1, l=0 state only
ground_state_data = data[(data['l'] == 0) & (data['n'] == 1)]

# Create a new figure with a reasonable size
plt.figure(figsize=(10, 6))

# Get unique iteration values (assuming they start from 0 and go to 19)
iterations = ground_state_data['iters'].unique()

# Create a colormap for the iterations
colors = plt.cm.viridis(np.linspace(0, 1, len(iterations)))

# Plot each iteration with a different color
for i, iter_num in enumerate(iterations):
    iter_data = ground_state_data[ground_state_data['iters'] == iter_num]
    # Sort by radius to ensure proper line plotting
    iter_data = iter_data.sort_values('r')
    plt.plot(iter_data['r'], iter_data['prob_density'], 
             label=f'Iteration {iter_num}', 
             color=colors[i], 
             linewidth=1.5)

# Add labels and title
plt.xlabel('Radial Distance (r)', fontsize=12)
plt.ylabel('Radial Probability Density', fontsize=12)
plt.title('Hydrogen-like Atom Radial Probability Density (n=1, l=0, Z=3)\nConvergence Over Iterations', fontsize=14)

# Add legend with smaller font and multiple columns if many iterations
plt.legend(fontsize=8, ncol=2)

# Set reasonable x-axis limit (adjust based on your data)
plt.xlim(0, 2)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Save the figure
plt.savefig('Plots/B3_4', dpi=300)
# plt.show()









import numpy as np
from scipy import integrate



# Load the radial probability density data
data = pd.read_csv('Hartree_Wavefunctions.txt', sep=' ', 
                   names=['iters','l', 'n', 'r', 'P'])


data = data[data['iters'] == number_iters]

initial_state = data[(data['l'] == 1) & (data['n'] == 2)]  # 2p state
final_state = data[(data['l'] == 0) & (data['n'] == 2)]    # 2s state


# Extract unique r values
r_values = initial_state['r'].unique()

# Create arrays for the wave functions
P_initial = initial_state['P'].values
P_final = final_state['P'].values

# Make sure the arrays have the same length
if len(r_values) == len(P_initial) and len(r_values) == len(P_final):

    # Calculate the integrand: P_initial * r * P_final
    integrand = P_initial * r_values * P_final

    # Perform the integration
    R_ab = np.trapz(integrand, r_values)

    omega_ab = 0.06791
    decay_rate = 2* R_ab**2 * omega_ab**3 / 3 * 1.071*pow(10,10)
    time = np.round(1/decay_rate * 10**9,4)
    print(time, " nanoseconds") 









