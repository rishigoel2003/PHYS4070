import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

def animation_SE_modified(file_name):
    # Load the data: columns are g, x, time, |ψ|²
    data = np.loadtxt(file_name)

    # Extract unique g values and time steps
    unique_g = np.unique(data[:, 0])
    unique_times = np.unique(data[:, 2])

    # Organize data per g and per time
    g_frames = {g: [] for g in unique_g}
    for g in unique_g:
        g_data = data[data[:, 0] == g]
        for t in unique_times:
            frame_data = g_data[g_data[:, 2] == t]
            frame_data = frame_data[frame_data[:, 1].argsort()]  # sort by x
            g_frames[g].append(frame_data)

    # Plot setup
    n_plots = len(unique_g)
    fig, axes = plt.subplots(n_plots, 1, figsize=(8, 3 * n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]  # make iterable

    lines = []
    texts = []

    x_min, x_max = np.min(data[:, 1]), np.max(data[:, 1])
    y_max = np.max(data[:, 3]) * 1.1

    for ax, g in zip(axes, unique_g):
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, y_max)
        ax.set_ylabel(r'$|\psi(x)|^2$')
        ax.set_title(f'g = {g}')
        ax.grid(True)
        line, = ax.plot([], [], lw=2)
        text = ax.text(0.85, 0.85, '', transform=ax.transAxes)
        lines.append(line)
        texts.append(text)

    axes[-1].set_xlabel('Position (x)')
    fig.tight_layout(pad=2.0)

    def animate(i):
        for idx, g in enumerate(unique_g):
            frame = g_frames[g][i]
            x = frame[:, 1]
            psi_sq = frame[:, 3]
            lines[idx].set_data(x, psi_sq)
            texts[idx].set_text(f'Time: {frame[0, 2]:.2f}')
        return lines + texts

    anim = FuncAnimation(fig, animate, frames=len(unique_times), interval=100, blit=True)

    # Save GIF
    gif_path = os.path.join(os.getcwd(), 'Plots/nlse_animation_plane_wave.gif')
    anim.save(gif_path, writer=PillowWriter(fps=10))

    # Save final frame as PNG
    # animate(-1)  # Draw final frame
    # png_path = os.path.join(os.getcwd(), 'Plots/nlse_final_frame_plane_wave.png')
    # fig.savefig(png_path)

    plt.close(fig)  # Close the figure to avoid duplicate display

# Call the function
animation_SE_modified("Part1/nlse_evolution_plane_wave.dat")
