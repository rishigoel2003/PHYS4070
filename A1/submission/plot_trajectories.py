import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import numpy as np




def plot_trajectories(filename):
    # Read the data from the text file
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

    
    # Convert lists to numpy arrays so we can plot nicely
    moon_x = np.array(moon_x)
    moon_y = np.array(moon_y)
    proj_x = np.array(proj_x)
    proj_y = np.array(proj_y)
    distances_moon_proj = np.array(distances_moon_proj)

    
    # Create figure for animation, setting so its a square frame with equal ticks
    fig_anim, ax_anim = plt.subplots(figsize=(8, 8))
    ax_anim.set_xlim(-20, 20)
    ax_anim.set_ylim(-20, 20)
    
    # Create line objects one for the moon and one for the projectile
    moon_line, = ax_anim.plot([], [], lw=2, color='b', label='Moon')
    proj_line, = ax_anim.plot([], [], lw=2, color='r', label='Projectile')
    
    # Add planet as a static circle (radius 1)
    planet = patches.Circle((0, 0), 1, linewidth=2, edgecolor='black', facecolor='none')
    ax_anim.add_patch(planet)
    
    # Add moon as a circle at the start point of the moon (radius 0.25)
    moon = patches.Circle((moon_x[0], moon_y[0]), 0.25, linewidth=2, edgecolor='blue', facecolor='none')
    ax_anim.add_patch(moon)
    
    #legend
    ax_anim.legend()
    #title
    ax_anim.set_title('Trajectory Animation')
    
    #setting the line to plot as the list of datapoints of the moon and projectiles up to a certain frame
    def animate(frame):
        moon_line.set_xdata(moon_x[:frame])
        moon_line.set_ydata(moon_y[:frame])
        proj_line.set_xdata(proj_x[:frame])
        proj_line.set_ydata(proj_y[:frame])
        #setting center of the moon to be the last datapoint
        moon.center = (moon_x[frame-1], moon_y[frame-1])
        return moon_line, proj_line, moon
    
    #creating animation using above function
    anim = FuncAnimation(fig_anim, animate, frames=len(moon_x), interval=20, blit=True)
    
    # Create static final trajectory plot 
    fig_final, ax_final = plt.subplots(figsize=(8, 8))
    ax_final.set_xlim(-20, 20)
    ax_final.set_ylim(-20, 20)
    
    #label the paths
    ax_final.plot(moon_x, moon_y, 'b-', linewidth=2, label='Moon Trajectory')
    ax_final.plot(proj_x, proj_y, 'r-', linewidth=2, label='Projectile Trajectory')
    
    # Add planet and final moon position as well as the final projectile
    ax_final.add_patch(patches.Circle((0, 0), 1, linewidth=2, edgecolor='black', facecolor='none'))
    ax_final.add_patch(patches.Circle((moon_x[-1], moon_y[-1]), 0.25, linewidth=2, edgecolor='blue', facecolor='none'))

    ax_final.plot(moon_x[-1], moon_y[-1], 'bo', markersize=8, label='Moon Final Position')
    ax_final.plot(proj_x[-1], proj_y[-1], 'ro', markersize=8, label='Projectile Final Position')
    
    #legend,title and grid lines
    ax_final.legend()
    ax_final.set_title('Complete Trajectories')
    ax_final.grid(True)




    ### === DISTANCE v TIME PLOT === ### used in later questions more so but good to plot for all the questions
    fig_dist = plt.figure(figsize=(8, 6))
    ax_dist = fig_dist.add_subplot(111)

    #plotting data
    ax_dist.plot(time, distances_moon_proj, 'g-', linewidth=2, label='Moon-Projectile Distance')

    # Add labels and title
    ax_dist.set_xlabel("Time")
    ax_dist.set_ylabel("Distance between Moon and Projectile")
    ax_dist.set_title("Moon-Projectile Distance vs. Time")
    ax_dist.grid(True)
    ax_dist.legend()


    ### === MOON PERSPECTIVE PLOT === ###
    # Create figure for final trajectory plot from moons perspective
    fig_final = plt.figure(figsize=(8, 8))
    ax_final = plt.axes(xlim=(-1, 1), ylim=(-1, 1))

    # Plot complete trajectories from moon's perspective
    ax_final.plot(proj_x-moon_x, proj_y-moon_y, 'r-', label='Projectile Trajectory - Moon Perspective', linewidth=2)

    # Add legend and title for final plot
    ax_final.legend()
    ax_final.set_title('Complete Trajectories moon perspective')
    ax_final.grid(True)




    
    plt.show()  # Show all figures at the same time
