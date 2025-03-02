import matplotlib.pyplot as plt
import pandas as pd

# Read the timing results
data = pd.read_csv('timing_results.txt')

# Create the plot
plt.figure(figsize=(10, 6))

# Plot the execution time vs number of threads
plt.plot(data['Threads'], data['Time'], 'o-', linewidth=2, markersize=8)
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (seconds)')
plt.title('Ising Model Simulation: Execution Time vs Number of Threads')
plt.grid(True)

# Add log2 scale for x-axis since we're using powers of 2
plt.xscale('log', base=2)
plt.xticks(data['Threads'], data['Threads'])

plt.show()
