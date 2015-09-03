__author__ = 'martin'

import matplotlib.pyplot as plt


def get_raw_peaks(path):
    peaks_file = open(path)
    peaks_string = peaks_file.read()
    peaks_file.close()
    peak_string_list = peaks_string.splitlines()
    peak_coordinates = [[float(y) for y in x.split("\t")[1:-1]] for x in peak_string_list]
    if [] in peak_coordinates:
        peak_coordinates.remove([])
    return peak_coordinates


def main():
    cosy_peak_points = get_raw_peaks('cosy peaks indanone.txt')
    print(cosy_peak_points)
    plot_graphs(cosy_peak_points)


def plot_graphs(coordinates):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title("X vs Y AVG",fontsize=14)
    ax.set_xlabel("XAVG",fontsize=12)
    ax.set_ylabel("YAVG",fontsize=12)
    ax.grid(True,linestyle='-',color='0.75')
    x_coords = [x[0] for x in coordinates]
    y_coords = [x[1] for x in coordinates]
    ax.scatter(x_coords,y_coords,s=20, marker='o')
    plt.show()

if __name__ == "__main__":
    main()
