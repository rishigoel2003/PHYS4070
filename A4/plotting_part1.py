import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Read the data file
data = np.loadtxt('nlse_evolution.dat')

# Group data by time steps
unique_times = np.unique(data[:, 1])
frames = []

for t in unique_times:
    # Get all rows with this time value
    frame_data = data[data[:, 1] == t]
    # Sort by x position to ensure proper order
    frame_data = frame_data[frame_data[:, 0].argsort()]
    frames.append(frame_data)

# Create the figure and axes
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
fig.tight_layout(pad=3.0)

# Initialize the plots
line_abs, = ax1.plot([], [], 'b-', lw=2, label=r'$|\psi(x)|$')
line_real, = ax2.plot([], [], 'r-', lw=2, label=r'Re$[\psi(x)]$')
line_imag, = ax2.plot([], [], 'g-', lw=2, label=r'Im$[\psi(x)]$')

# Configure the axes
ax1.set_title('NLSE: Bright Soliton Evolution')
ax1.set_ylabel('Amplitude')
ax1.grid(True)
ax1.legend()

ax2.set_xlabel('Position (x)')
ax2.set_ylabel('Re/Im Parts')
ax2.grid(True)
ax2.legend()

# Set the axis limits
x_min, x_max = np.min(data[:, 0]), np.max(data[:, 0])
abs_max = np.max(data[:, 2]) * 1.1
re_im_min = np.min(data[:, 3:5]) * 1.1
re_im_max = np.max(data[:, 3:5]) * 1.1

ax1.set_xlim(x_min, x_max)
ax1.set_ylim(0, abs_max)
ax2.set_xlim(x_min, x_max)
ax2.set_ylim(re_im_min, re_im_max)

# Add time indicator
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

def animate(i):
    frame = frames[i]
    x_data = frame[:, 0]
    
    line_abs.set_data(x_data, frame[:, 2])  # |ψ|
    line_real.set_data(x_data, frame[:, 3])  # Re[ψ]
    line_imag.set_data(x_data, frame[:, 4])  # Im[ψ]
    
    time_text.set_text(f'Time: {frame[0, 1]:.2f}')
    
    return line_abs, line_real, line_imag, time_text

# Create the animation
anim = FuncAnimation(
    fig, animate, frames=len(frames),
    interval=200, blit=True
)

plt.show()