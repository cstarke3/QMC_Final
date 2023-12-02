import matplotlib.pyplot as plt
import numpy as np

def mean(lst, pct_val = 0.3):
    
    # mean_val = mean_last_r(ys, int(len(ys)*pct_val)) # get the mean using the last 10% of the data
    r = int(len(lst)*pct_val)
    mean_val = sum(lst[-r:]) / r
    return mean_val

def mean_stddev(lst, pct_val = 0.5):
    start_index = int(len(lst) * (1 - pct_val))

    # Slice the array to get the last X percent
    last_percent_data = lst[start_index:]

    # Compute mean and standard deviation
    mean = np.mean(last_percent_data)
    stddev = np.std(last_percent_data)
    return mean, stddev
    
def plot_data(ys, n_vals, hist_data, bins, title="QMC Simulation"):
    """ Plots the data with an inset histogram. """

    mean_val, stddev = mean_stddev(ys)

    xs = range(len(ys))

    # Create a new figure for the main plot
    plt.figure()


    # Plot x vs y
    plt.plot(xs, ys, label="E_ref")
    plt.xlabel('Time Step')
    plt.ylabel('E_ref')
    plt.title(title, fontsize=16)
    plt.hlines(mean_val, 0, len(ys), colors='r', linestyles='solid', label='E_0')
    text_x_position = 0.95 * len(xs)  # X position near the right edge of the plot
    text_y_position = mean_val + 0.03  # Y position for the text just above the line
    plt.text(text_x_position, text_y_position, 
             f"E_0: {mean_val:.3f} +/- {stddev:.3f}", 
             fontsize=14, color='black', horizontalalignment='right',
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    # Plot the replica count on a secondary axis (in the background)
    # ax_twin = plt.twinx()
    # ax_twin.plot(xs, n_vals, label="Replica Count", color='green', alpha=0.3)
    plt.legend(loc='lower right')
    # Create an inset plot for the histogram
    # The arguments are [left, bottom, width, height] in figure coordinate
    plot_inset = False
    if plot_inset:
        ax_inset = plt.axes([0.45, 0.2, 0.4, 0.4])  # Modify these values to adjust the position and size of the inset
        ax_inset.hist(hist_data, bins=bins, color='green', alpha=0.5)
        ax_inset.set_title("Ground State Wave Function")
        ax_inset.set_xlabel('Position')
        ax_inset.set_ylabel('Count')

    # Show the plot
    plt.show()
    

def plot_energy_vs_alpha(alpha_xs, energy_ys, title="Energy vs Alpha"):
    """ Plots the data with an inset histogram. """

    mean_val, stddev = mean_stddev(energy_ys)

    # Create a new figure for the main plot
    plt.figure()

    xmin = 0.0
    xmax = 1.0
    # Plot x vs y
    plt.plot(alpha_xs, energy_ys, label="E_0_mean")
    plt.xlabel('alpha vals')
    # plt.xlim(alpha_min, alpha_max)
    plt.ylabel('E_0_mean')
    plt.title(title, fontsize=16)
    plt.hlines(mean_val, xmin, xmax, colors='r', linestyles='solid', label='E_0')
    text_x_position = 0.95 * len(alpha_xs)  # X position near the right edge of the plot
    text_y_position = mean_val + 0.03  # Y position for the text just above the line
    plt.text(text_x_position, text_y_position, 
             f"E_0: {mean_val:.3f} +/- {stddev:.3f}", 
             fontsize=14, color='black', horizontalalignment='right',
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    # Plot the replica count on a secondary axis (in the background)
    # ax_twin = plt.twinx()
    # ax_twin.plot(xs, n_vals, label="Replica Count", color='green', alpha=0.3)
    plt.legend(loc='lower right')
 
    # Show the plot
    plt.show()

    
def plot_histogram(data, bins, title = "QMC Simulation"):
    """ Plots the data. """
    
    # Create a new figure
    plt.figure()

    # Plot x vs y
    plt.hist(data, bins=bins, label="positions")

    # Add title and labels
    plt.title(title, fontsize=16)
    plt.xlabel('Centroid Positions')
    plt.ylabel('Count')

    # Add a legend
    # plt.legend()

    # Show the plot
    plt.show()
    