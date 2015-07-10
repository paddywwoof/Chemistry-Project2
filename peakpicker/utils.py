import matplotlib.pyplot as plt
import os

def readfile(path):
    current_path = os.path.dirname(__file__) + "/"
    openfile = open(current_path+path, "r")
    file_string = openfile.read()
    openfile.close()
    return file_string


def plot_graph(peaks):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.grid(True,linestyle='-',color='0.75')
    x_coords = [x[0] for x in peaks]
    y_coords = [x[1] for x in peaks]
    ax.scatter(x_coords,y_coords,s=20, marker='o')
    plt.show()
