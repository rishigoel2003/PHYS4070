import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

def animation_SE_modified(file_name):
    # Load the data: columns are u, x, time, |ψ|²
    data = np.loadtxt(file_name)

    # Extract unique u values and time steps
    unique_u = np.unique(data[:, 0])
    unique_times = np.unique(data[:, 2])

    # Organize data per u and per time
    u_frames = {u: [] for u in unique_u}
    for u in unique_u:
        u_data = data[data[:, 0] == u]
        for t in unique_times:
            frame_data = u_data[u_data[:, 2] == t]
            frame_data = frame_data[frame_data[:, 1].argsort()]  # sort by x
            u_frames[u].append(frame_data)

    # Plot setup
    n_plots = len(unique_u)
    fig, axes = plt.subplots(n_plots, 1, figsize=(8, 3 * n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]  # make iterable

    lines = []
    texts = []

    x_min, x_max = np.min(data[:, 1]), np.max(data[:, 1])
    y_max = np.max(data[:, 3]) * 1.1

    for ax, u in zip(axes, unique_u):
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, y_max)
        ax.set_ylabel(r'$|\psi(x)|^2$')
        ax.set_title(f'u = {u}')
        ax.grid(True)
        line, = ax.plot([], [], lw=2)
        text = ax.text(0.85, 0.85, '', transform=ax.transAxes)
        lines.append(line)
        texts.append(text)

    axes[-1].set_xlabel('Position (x)')
    fig.tight_layout(pad=2.0)

    def animate(i):
        for idx, u in enumerate(unique_u):
            frame = u_frames[u][i]
            x = frame[:, 1]
            psi_sq = frame[:, 3]
            lines[idx].set_data(x, psi_sq)
            texts[idx].set_text(f'Time: {frame[0, 2]:.2f}')
        return lines + texts

    anim = FuncAnimation(fig, animate, frames=len(unique_times), interval=100, blit=True)

    # Save GIF
    gif_path = os.path.join(os.getcwd(), 'Plots/nlse_animation_wave_packet.gif')
    anim.save(gif_path, writer=PillowWriter(fps=10))

    # Save final frame as PNG
    # animate(-1)  # Draw final frame
    # png_path = os.path.join(os.getcwd(), 'Plots/nlse_final_frame_wave_packet.png')
    # fig.savefig(png_path)

    plt.close(fig)  # Close the figure to avoid duplicate display

# Call the function
animation_SE_modified("Part1/nlse_evolution_wave_packet.dat")










import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt("Part1/nlse_evolution_wave_packet_peaks.dat")

# Extract u values
u_values = np.unique(data[:, 0])

plt.figure(figsize=(8, 5))

# Loop through u and plot corresponding t vs x_peak
for u in u_values:
    mask = data[:, 0] == u
    t_vals = data[mask, 1]
    x_peak_vals = data[mask, 2]
    plt.plot(t_vals, x_peak_vals, label=f'u = {u:.2f}')

plt.xlabel('Time')
plt.ylabel('Peak position $x_\\mathrm{peak}$')
plt.title('Wave Packet Peak Evolution for Different $u$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Plots/nlse_wave_packet_peaks.png")
# plt.show()
