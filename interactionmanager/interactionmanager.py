__author__ = 'martin'
import math
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt


class Interaction:
    def __init__(self, interaction_name, x_axis, repulsive_amplitude, repulsive_time_constant, attractive_amplitude,
                 attractive_time_constant, depth, power):
        self.interaction_name = interaction_name
        distance_function = [0 for i in x_axis]

        # To make equation easier to read
        aad = attractive_amplitude * depth
        rad = repulsive_amplitude * depth
        atc = attractive_time_constant
        rtc = repulsive_time_constant

        for xi in x_axis:
            distance_function[xi] = ((-aad * math.exp(-xi/atc)) + aad + (rad * math.exp(-xi/rtc)) - rad)**power
            if distance_function[xi] > distance_function[xi-1] and distance_function[xi-1] < distance_function[xi-2]:
                print(self.interaction_name+", Minimum at: ", xi-1)
        self.distance_function = distance_function

    def plot_graph(self, x_axis):
        number_points = 400
        plt.plot(x_axis[:number_points], self.distance_function[:number_points])
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.title('Two Dimensional')
        print("Plotting: " + self.interaction_name)

    def get_response_value(self,distance):
        return self.distance_function[distance]

class DefaultInteraction(Interaction):
    def __init__(self, x_axis, repulsive_amplitude, repulsive_time_constant, depth):
        self.interaction_name = "Default"
        self.depth = depth
        self.distance_function = [0 for i in x_axis]
        for xi in x_axis:
            self.distance_function[xi] = self.depth * (repulsive_amplitude * math.exp(-xi/repulsive_time_constant))

    def get_response_value(self, distance):
        return self.distance_function[distance]

class InteractionManager:
    """
    Managers the different types of interaction (e.g. HSQC, COSY etc, and their response values)
    """
    def __init__(self, axis_width, interaction_matrix, atom_types, shift_data):
        """
        Args:
            axis_width: Number of points to generate function mapping for

        Returns:
            None
        """
        self.x_axis = range(axis_width)
        self.interaction_map = {}
        self.interaction_matrix = interaction_matrix
        self.atom_types = atom_types
        self.shift_data = shift_data
        self.number_carbon_signals = self.atom_types.count("C")
        self.number_hydrogen_signals = self.atom_types.count("H")
        self.number_signals = self.number_carbon_signals + self.number_hydrogen_signals
        self.bonds = self.get_bonds()

    def add_default_interaction(self,index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth):
        """
        Adds default repulsive interaction type
        """
        new_interaction = DefaultInteraction(self.x_axis, repulsive_amplitude, repulsive_time_constant, depth)
        self.interaction_map[index] = new_interaction

    def add_new_interaction(self, index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth,
                            attractive_amplitude, attractive_time_constant, power):
        """
        Adds a new interaction to the interaction manager

        Args:
            interaction_name: Name of the type of interaction (E.g. COSY, HSQC etc)
            repulsive amplitude : Repulsive Amplitude
            repulsive_time_constant : Repulsive Time Constant
            depth : Depth
            attractive_amplitude: Attractive Amplitude
            power : Power

        Returns:
            None
        """
        new_interaction = Interaction(interaction_name, self.x_axis, repulsive_amplitude, repulsive_time_constant, attractive_amplitude,
                                      attractive_time_constant, depth, power)
        self.interaction_map[index] = new_interaction

    def calculate_magnitude(self, v1,v2):
        distance = sqrt( (v1[0]-v2[0]) ** 2 + (v1[1]-v2[1]) ** 2 + (v1[2]-v2[2]) ** 2)

        return min(1000,round(distance))

    def interaction_response(self, i, j, v1, v2):
        """
        Args:
            i  : ith Atom index in the array
            j  : jth Atom index in the array
            v1 : Atom 1
            v2 : Atom 2
        Returns:
            interaction response
        """
        try:
            distance = self.calculate_magnitude(v1, v2)
            interaction_type = self.interaction_matrix[j][i]
            interaction = self.interaction_map[interaction_type]
            return interaction.get_response_value(distance)
        except IndexError:
            return 0
        except KeyError:
            return 0

    def shape_coordinates(self,atom_coordinates):
        """
        When array is improperly shaped, reshapes the array
        """
        return atom_coordinates.reshape((self.number_signals, 3))

    def plot_all_interactions(self):

        for interaction in self.interaction_map.values():
            interaction.plot_graph(self.x_axis)
        #pyplot.ion()
        plt.show()

    def get_initial_coordinates(self):
        atom_coordinates = np.random.uniform(99, 101, (self.number_signals, 3))
        return atom_coordinates

    def between(self, x, interval):
        if interval[0] <= x < interval[1]:
            return True
        else:
            return False

    def get_bonds(self):
        singlebond = [0,100]
        doublebond = [100,160]
        oxygen = [160,222]

        bonds = set()
        for i in range(len(self.interaction_matrix)):
            for j in range(len(self.interaction_matrix)):
                if self.interaction_matrix[i][j] in [2, 4]:
                    if self.between(self.shift_data[i], singlebond) and self.between(self.shift_data[j], singlebond):
                        btype = 1
                    if self.between(self.shift_data[i], doublebond) and self.between(self.shift_data[j], doublebond):
                        btype = 2
                    else:
                        btype = 1
                    new_bond = [i+1, j+1]
                    new_bond.sort()
                    bonds.add(tuple(new_bond+[btype]))
        bonds = list(bonds)
        #Add Oxygens
        for shift_index,shift_value in enumerate(self.shift_data):
            if self.between(shift_value, oxygen):
                self.atom_types.append("O")
                bonds.append([shift_index+1,len(self.atom_types),2])
        return bonds

    def print_matrix(self):
        print("     ",*[str(x+1)+" "*(3-len(str(x+1))) for x in range(self.number_signals)])
        for x in range(self.number_signals):
            print(str(x+1)+" "*(2-len(str(x+1))), self.atom_types[x], str(self.interaction_matrix[x]).replace(" ", "   "), self.shift_data[x])


