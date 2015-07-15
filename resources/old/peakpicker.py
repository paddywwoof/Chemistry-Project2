import numpy as np
import matplotlib.pyplot as plt

from old.utils import round_to_nearest, merge_duplicates, remove_asymmetric


class PeakManager:
    def __init__(self, path):
        self.peak_points = self.read_peaks_from_file(path)
        self.diagonals = self.get_diagonal_peaks()

    def read_peaks_from_file(self, path):
        """
        Args:
            path : path of peaks text file
        Returns:
            peak_coordinates : raw nmr peak data
        """
        peaks_file = open(path)
        peaks_string = peaks_file.read()
        peaks_file.close()
        peak_string_list = peaks_string.splitlines()
        peak_points = [[float(y) for y in x.split("\t")[1:-1]] for x in peak_string_list]
        if [] in peak_points:
            peak_points.remove([])
        return peak_points

    def get_diagonal_peaks(self):
        """
        Finds all points that approximate diagonals (when the coordinates are almost the same)
        and generates a list of these points.
        Similar diagonals are then merged together

        Args:
            peak_points : raw nmr peak data
        Returns:
            grouped_diagonals : list of rounded diagonal peak points
        """
        error = 0.05  # amount that x can vary from y e.g [3.65,3.55] fails, [3.65,3.63] passes
        diagonals_list = []
        for peak in self.peak_points:
            if abs(peak[0]-peak[1]) < error and peak[2]>0.3:
                diag_value = round(0.5*peak[0]+0.5*peak[1], 2)
                diagonals_list.append([diag_value, diag_value, peak[2]])

        current_group = []
        grouped_diagonals = []
        for diag in diagonals_list:
            if not current_group:
                current_group.append(diag)
            else:
                if abs(current_group[-1][0] - diag[0]) < 0.075:
                    current_group.append(diag)
                else:
                    grouped_diagonals.append(self.merge_group(current_group))
                    current_group = [diag]
        grouped_diagonals.append(self.merge_group(current_group))
        return grouped_diagonals

    def merge_group(self, group):
        weight = sum([x[2] for x in group])
        weighted_x = round( sum([x[0]*x[2]/weight for x in group]) , 2)
        return [weighted_x, weighted_x, round(weight,6)]

    def refine_peak_points(self, peak_points, diagonals_list):
        """
        Pushes all peaks to nearest hydrogen coordinates from diagonals
        If it varies in both coordinates by more than 0.1, invalid coordinate

        Args:
            peak_points: list of all unrefined peak points
            diagonals_list: list of diagonal peak points
        Returns:
            new_peak_points: refined peak points

        """
        error = 0.1
        new_peak_points = []
        for peak_index,peak in enumerate(peak_points):
            new_peak = [round_to_nearest(diagonals_list, peak[0]), round_to_nearest(diagonals_list, peak[1]), peak[2]]
            if abs(peak[0]-new_peak[0])<0.1 and abs(peak[1]-new_peak[1]) < error:
                new_peak_points.append(new_peak)

        new_peak_points = merge_duplicates(new_peak_points)
        new_peak_points.sort(key=lambda x: x[0])
        return new_peak_points

    def main(self):
        diagonals = self.diagonals
        diagonal_values = [x[0] for x in diagonals]

        rounded_peak_points = self.refine_peak_points(self.peak_points, [x[0] for x in diagonals])
        rounded_peak_points = [x for x in rounded_peak_points if [x[1], x[0]] in [y[0:2] for y in rounded_peak_points] ]
        print(rounded_peak_points)
        plot_graphs(rounded_peak_points)
        cosy_interactions = [[diagonal_values.index(x[0]), diagonal_values.index(x[1])] for x in rounded_peak_points]
        cosy_interactions = remove_asymmetric(cosy_interactions)
        print(cosy_interactions)
        interaction_matrix = np.zeros((len(diagonal_values), len(diagonal_values)),dtype=np.int)
        for i in cosy_interactions:
            if i[0]!=i[1]:
                interaction_matrix[i[0]][i[1]] = 1
                interaction_matrix[i[1]][i[0]] = 1
        print(interaction_matrix)


def plot_graphs(peaks):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.grid(True,linestyle='-',color='0.75')
    x_coords = [x[0] for x in peaks]
    y_coords = [x[1] for x in peaks]
    ax.scatter(x_coords,y_coords,s=20, marker='o')
    plt.show()

if __name__ == "__main__":
    PeakManager('cosy peaks indanone.txt').main()


