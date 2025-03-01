import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import numpy as np




def plot_trajectories(filename):
    # Read the data from the file
    time, moon_x, moon_y, proj_x, proj_y, distances_moon_proj = [], [], [], [], [], []
    
    with open(filename, "r") as file:
        for line in file:
            data = line.split()
            time.append(float(data[0]))
            moon_x.append(float(data[1]))
            moon_y.append(float(data[2]))
            proj_x.append(float(data[3]))
            proj_y.append(float(data[4]))
            distances_moon_proj.append(float(data[5]))

    
    # Convert lists to numpy arrays
    moon_x = np.array(moon_x)
    moon_y = np.array(moon_y)
    proj_x = np.array(proj_x)
    proj_y = np.array(proj_y)
    distances_moon_proj = np.array(distances_moon_proj)

    
    # Create figure for animation
    fig_anim, ax_anim = plt.subplots(figsize=(8, 8))
    ax_anim.set_xlim(-20, 20)
    ax_anim.set_ylim(-20, 20)
    
    # Create line objects
    moon_line, = ax_anim.plot([], [], lw=2, color='b', label='Moon')
    proj_line, = ax_anim.plot([], [], lw=2, color='r', label='Projectile')
    
    # Add Earth
    earth = patches.Circle((0, 0), 1, linewidth=2, edgecolor='black', facecolor='none')
    ax_anim.add_patch(earth)
    
    # Add moon
    moon = patches.Circle((moon_x[0], moon_y[0]), 0.25, linewidth=2, edgecolor='blue', facecolor='none')
    ax_anim.add_patch(moon)
    
    ax_anim.legend()
    ax_anim.set_title('Trajectory Animation')
    
    def animate(frame):
        moon_line.set_xdata(moon_x[:frame])
        moon_line.set_ydata(moon_y[:frame])
        proj_line.set_xdata(proj_x[:frame])
        proj_line.set_ydata(proj_y[:frame])
        moon.center = (moon_x[frame-1], moon_y[frame-1])
        return moon_line, proj_line, moon
    
    anim = FuncAnimation(fig_anim, animate, frames=len(moon_x), interval=20, blit=True)
    
    # Create final trajectory plot
    fig_final, ax_final = plt.subplots(figsize=(8, 8))
    ax_final.set_xlim(-20, 20)
    ax_final.set_ylim(-20, 20)
    
    ax_final.plot(moon_x, moon_y, 'b-', linewidth=2, label='Moon Trajectory')
    ax_final.plot(proj_x, proj_y, 'r-', linewidth=2, label='Projectile Trajectory')
    
    # Add Earth and final moon position
    ax_final.add_patch(patches.Circle((0, 0), 1, linewidth=2, edgecolor='black', facecolor='none'))
    ax_final.add_patch(patches.Circle((moon_x[-1], moon_y[-1]), 0.25, linewidth=2, edgecolor='blue', facecolor='none'))
    
    ax_final.plot(moon_x[-1], moon_y[-1], 'bo', markersize=8, label='Moon Final Position')
    ax_final.plot(proj_x[-1], proj_y[-1], 'ro', markersize=8, label='Projectile Final Position')
    
    ax_final.legend()
    ax_final.set_title('Complete Trajectories')
    ax_final.grid(True)




    ### === DISTANCE vs. TIME PLOT === ###
    fig_dist = plt.figure(figsize=(8, 6))
    ax_dist = fig_dist.add_subplot(111)

    ax_dist.plot(time, distances_moon_proj, 'g-', linewidth=2, label='Moon-Projectile Distance')

    # Add labels and title
    ax_dist.set_xlabel("Time")
    ax_dist.set_ylabel("Distance between Moon and Projectile")
    ax_dist.set_title("Moon-Projectile Distance vs. Time")
    ax_dist.grid(True)
    ax_dist.legend()


    
    # Create figure for final trajectory plot from moons perspective
    fig_final = plt.figure(figsize=(8, 8))
    ax_final = plt.axes(xlim=(-1, 1), ylim=(-1, 1))

    # Plot complete trajectories
    # ax_final.plot(moon_x, moon_y, 'b-', label='Moon Trajectory', linewidth=2)
    ax_final.plot(proj_x-moon_x, proj_y-moon_y, 'r-', label='Projectile Trajectory - Moon Perspective', linewidth=2)


    # ax_final.plot(proj_x[-1]-moon_x[-1], proj_y[-1]-moon_y[-1], 'ro', markersize=8, label='Projectile Final Position')

    # Add legend and title for final plot
    ax_final.legend()
    ax_final.set_title('Complete Trajectories moon perspective')
    ax_final.grid(True)




    
    plt.show()  # Show both figures at the same time
