import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

def animation_SE_modified(file_name):
    # Load the data: columns are theta, x, time, |ψ|²
    data = np.loadtxt(file_name)

    # Extract unique theta values and time steps
    unique_theta = np.unique(data[:, 0])
    unique_times = np.unique(data[:, 2])

    # Organize data per theta and per time
    theta_frames = {theta: [] for theta in unique_theta}
    for theta in unique_theta:
        theta_data = data[data[:, 0] == theta]
        for t in unique_times:
            frame_data = theta_data[theta_data[:, 2] == t]
            frame_data = frame_data[frame_data[:, 1].argsort()]  # sort by x
            theta_frames[theta].append(frame_data)

    # Plot setup
    n_plots = len(unique_theta)
    fig, axes = plt.subplots(n_plots, 1, figsize=(8, 3 * n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]  # make iterable

    lines = []
    texts = []

    x_min, x_max = np.min(data[:, 1]), np.max(data[:, 1])
    y_max = np.max(data[:, 3]) * 1.1

    for ax, theta in zip(axes, unique_theta):
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, y_max)
        ax.set_ylabel(r'$|\psi(x)|^2$')
        ax.set_title(f'theta = {theta}')
        ax.grid(True)
        line, = ax.plot([], [], lw=2)
        text = ax.text(0.85, 0.85, '', transform=ax.transAxes)
        lines.append(line)
        texts.append(text)

    axes[-1].set_xlabel('Position (x)')
    fig.tight_layout(pad=2.0)

    def animate(i):
        for idx, theta in enumerate(unique_theta):
            frame = theta_frames[theta][i]
            x = frame[:, 1]
            psi_sq = frame[:, 3]
            lines[idx].set_data(x, psi_sq)
            texts[idx].set_text(f'Time: {frame[0, 2]:.2f}')
        return lines + texts

    anim = FuncAnimation(fig, animate, frames=len(unique_times), interval=100, blit=True)

    # Save GIF
    gif_path = os.path.join(os.getcwd(), 'Plots/nlse_animation_solitons.gif')
    anim.save(gif_path, writer=PillowWriter(fps=10))

    # Save final frame as PNG
    # animate(-1)  # Draw final frame
    # png_path = os.path.join(os.getcwd(), 'Plots/nlse_final_frame_solitons.png')
    # fig.savefig(png_path)

    plt.close(fig)  # Close the figure to avoid duplicate display

# Call the function
animation_SE_modified("Part1/nlse_evolution_solitons.dat")








import matplotlib.pyplot as plt
import numpy as np
import os

# Load data
data = np.loadtxt("Part1/nlse_evolution_solitons_peaks.dat")

# Unique theta values (rounded for stability)
theta_values = np.unique(np.round(data[:, 0], 8))  # Reduce floating-point noise

# Helper function for stable masking
def mask_for_theta(data, theta, tol=1e-6):
    return np.isclose(data[:, 0], theta, atol=tol)

# Create output directory
os.makedirs("Plots", exist_ok=True)

# === Plot 1: Subplots of Peak Positions ===
fig1, axes1 = plt.subplots(len(theta_values), 1, figsize=(8, 4 * len(theta_values)), sharex=True)
for i, theta in enumerate(theta_values):
    ax = axes1[i] if len(theta_values) > 1 else axes1
    mask = mask_for_theta(data, theta)
    t = data[mask, 1]
    x1 = data[mask, 2]
    x2 = data[mask, 3]
    ax.plot(t, x1, label="$x_1$", color='blue')
    ax.plot(t, x2, label="$x_2$", color='orange', linestyle='--')
    ax.set_ylabel("Position")
    ax.set_title(f"θ = {theta:.2f}")
    ax.grid(True)
    ax.legend()
axes1[-1].set_xlabel("Time")
plt.tight_layout()
plt.savefig("Plots/nlse_peaks_vs_time_subplots.png")

## === Plot: All Peak Distances on One Plot ===
fig, ax = plt.subplots(figsize=(10, 6))
for theta in theta_values:
    mask = mask_for_theta(data, theta)
    t = data[mask, 1]
    d = data[mask, 4]
    ax.plot(t, d, label=f"θ = {theta:.2f}")
ax.set_xlabel("Time")
ax.set_ylabel("Distance")
ax.set_title("Peak Distance vs Time for Different θ")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.savefig("Plots/nlse_peak_distance_vs_time.png")


# === Plot 3: Min Peak Distance vs θ ===
min_distances = []
for theta in theta_values:
    mask = mask_for_theta(data, theta)
    d = data[mask, 4]
    min_distances.append(np.min(d))

plt.figure(figsize=(8, 5))
plt.plot(theta_values, min_distances, 'o-', color='red')
plt.xlabel('θ')
plt.ylabel('Minimum Peak Distance')
plt.title('Minimum Peak Distance vs θ')
plt.grid(True)
plt.tight_layout()
plt.savefig("Plots/nlse_min_peak_distance_vs_theta.png")

