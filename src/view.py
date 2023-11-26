import matplotlib.pyplot as plt

def mean(lst, pct_val = 0.3):
    
    # mean_val = mean_last_r(ys, int(len(ys)*pct_val)) # get the mean using the last 10% of the data
    r = int(len(lst)*pct_val)
    mean_val = sum(lst[-r:]) / r
    return mean_val

def plot_data(ys, data, bins, title="QMC Simulation"):
    """ Plots the data with an inset histogram. """

    mean_val = mean(ys)
    xs = range(len(ys))

    # Create a new figure for the main plot
    plt.figure()

    # Plot x vs y for the main plot
    plt.plot(xs, ys, label="E_ref")
    plt.xlabel('Time Step')
    plt.ylabel('E_ref')
    plt.title(title, fontsize=16)
    plt.hlines(mean_val, 0, len(ys), colors='r', linestyles='solid', label='E_0')
    text_x_position = 0.95 * len(xs)  # X position near the right edge of the plot
    text_y_position = mean_val + 0.03  # Y position for the text just above the line
    plt.text(text_x_position, text_y_position, 
             f"Mean: {mean_val:.2f}", fontsize=14, color='black', horizontalalignment='right',
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    # Create an inset plot for the histogram
    # The arguments are [left, bottom, width, height] in figure coordinate
    ax_inset = plt.axes([0.45, 0.2, 0.4, 0.4])  # Modify these values to adjust the position and size of the inset
    ax_inset.hist(data, bins=bins, color='green', alpha=0.5)
    ax_inset.set_title("Ground State Wave Function")
    ax_inset.set_xlabel('Position')
    ax_inset.set_ylabel('Count')

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
    