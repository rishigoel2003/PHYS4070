import pandas as pd
import numpy as np
import matplotlib.pyplot as plt





# Load data
data = pd.read_csv('Greens_Energy.txt', sep=' ', names=['l', 'n', 'energy'])

# Filter for n = 2 only
n_2_data = data[data['n'] == 2]

# Theoretical/experimental energy values for 2s and 2p in atomic units
theoretical_energies = [-0.19814, -0.13023]

# Plot the data
plt.figure(figsize=(8, 5))

# Define colors for different l values
colors = {0: 'green', 1: 'blue'}

# Plot simulation values for different l values at n=2
unique_ls = n_2_data['l'].unique()
for l in unique_ls:
    l_data = n_2_data[n_2_data['l'] == l]
    plt.scatter(l_data['n'], l_data['energy'], color=colors[l], label=f'Simulation Energy (l={l})', marker='o')

# Plot theoretical values for 2s and 2p
plt.axhline(y=theoretical_energies[0], color='green', linestyle='--', label='Experimental 2s Energy')
plt.axhline(y=theoretical_energies[1], color='blue', linestyle='--', label='Experimental 2p Energy')

# Labels and title
plt.xlabel('Quantum number n')
plt.xticks([2])
plt.ylabel('Energy (au)')
plt.title('Simulation vs. Experimental Energy for 2s and 2p (Greens Approx)')
plt.legend()

# Show plot
plt.grid(True)
plt.savefig('Plots/B2_1', dpi=300)
# plt.show()













































# Load the radial probability density data
data = pd.read_csv('Greens_Probability_Density.txt', sep=' ', 
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
plt.title('Greens Potential - Radial Probability Densities (Z=3)', fontsize=14)

# Add legend
plt.legend(fontsize=10)

# Set reasonable x-axis limit (adjusted based on data)
plt.xlim(0, 15)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Save the figure
plt.savefig('Plots/B2_2', dpi=300)
# plt.show()












import numpy as np
from scipy import integrate



# Load the radial probability density data
data = pd.read_csv('Greens_Wavefunctions.txt', sep=' ', 
                   names=['l', 'n', 'r', 'P'])

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

#im doing this in python to check the result of the integration in c++ (it should give the same answer, i do the same 
#for b3 and b4)











delta = pd.read_csv('Greens_delta_E.txt', sep=' ', 
                   names=['l', 'de'])



# Load data
data = pd.read_csv('Greens_Energy.txt', sep=' ', names=['l', 'n', 'energy'])

# Extract relevant columns
ls = data['l']
ns = data['n']
energies = data['energy']

# Filter for n = 2 only
n_2_data = data[data['n'] == 2]

# Plot the data
plt.figure(figsize=(8, 5))

# Define colors for different l values
colors = {0: 'green', 1: 'blue'}


# Plot simulation values for n=2
unique_ls = n_2_data['l'].unique()
for l in unique_ls:
    l_data = n_2_data[n_2_data['l'] == l]

    delta_correction = delta.loc[delta['l'] == l, 'de'].values
    delta_correction = delta_correction[0]  # Extract scalar value
    # print((theoretical_energies[l]-l_data['energy']- delta_correction).values) (this was for me to see the difference and check if it works)

    plt.scatter(l_data['n'], l_data['energy'] + delta_correction, color = colors[l],
                label=f'simulation Energy (l={l})', marker='o')
    

# Plot experimental values for 2s and 2p
plt.axhline(y=theoretical_energies[0], color='green', linestyle='--', label='experimental 2s Energy')
plt.axhline(y=theoretical_energies[1], color='blue', linestyle='--', label='experimental 2p Energy')

# Labels and title
plt.xlabel('Quantum number n')
plt.ylabel('Energy (au)')
plt.title('simulation vs. experimental Energy with perturbation theory Correction')
plt.legend()

# Show plot
plt.grid(True)
plt.savefig('Plots/B2_3', dpi=300)
# plt.show()
