import matplotlib.pyplot as plt

def mean(lst, pct_val = 0.3):
    
    # mean_val = mean_last_r(ys, int(len(ys)*pct_val)) # get the mean using the last 10% of the data
    r = int(len(lst)*pct_val)
    mean_val = sum(lst[-r:]) / r
    return mean_val

def plot_data(xs,ys, title = "QMC Simulation"):
    """ Plots the data. """
    
    mean_val = mean(ys)
    # Create a new figure
    plt.figure()

    # Plot x vs y
    plt.plot(xs, ys, label="E_ref")

    # Add title and labels
    plt.title(title, fontsize=16)
    plt.xlabel('Time Step')
    plt.ylabel('E_ref')
    plt.hlines(mean_val, 0, len(ys), colors='r', linestyles='solid', label=f'E_0')
    text_x_position = 5  # X position where you want the text
    text_y_position = mean_val + 0.03  # Adjust text placement just above the line
    plt.text(text_x_position, text_y_position, f"Mean: {mean_val:.2f}", fontsize=14, color='black')

    # Add a legend
    plt.legend()

    # Show the plot
    plt.show()    