__author__ = 'martin'
import numpy as np
import math
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
            if distance_function[xi]>distance_function[xi-1] and distance_function[xi-1] < distance_function[xi-2]:
                print(self.interaction_name+", Minima at: ", xi)
        self.distance_function = distance_function

    def plot_graph(self, x_axis):
        number_points = 400
        plt.plot(x_axis[:number_points], self.distance_function[:number_points])
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.title('Two Dimensional')
        print("Plotting: " + self.interaction_name)


class DefaultInteraction(Interaction):
    def __init__(self, x_axis, repulsive_amplitude, repulsive_time_constant, depth):
        self.interaction_name = "Default"
        Default = [0 for i in x_axis]
        for xi in x_axis:
            Default[xi] = (repulsive_amplitude*depth * math.exp(-xi/repulsive_time_constant))
        self.distance_function = Default


class InteractionManager:
    """
    Managers the different types of interaction (e.g. HSQC, COSY etc, and their response values)
    """
    def __init__(self, get_interactions, axis_width):
        """
        Args:
            axis_width: Number of points to generate function mapping for

        Returns:
            None
        """

        self.x_axis = range(axis_width)
        self.interaction_map = {}
        self.interaction_matrix = get_interactions()
        self.atom_types = self.get_atom_types()
        self.number_carbon_signals = self.atom_types.count("C")
        self.number_hydrogen_signals= self.atom_types.count("H")
        self.number_signals = self.number_carbon_signals + self.number_hydrogen_signals


    def add_default_interaction(self,index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth):
        """
        Adds default repulsive interaction type
        """
        new_interaction = DefaultInteraction(self.x_axis, repulsive_amplitude, repulsive_time_constant,depth)
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

    def interaction_response(self, i, j, distance):
        """
        Args:
            i : ith Atom in the array
            j : jth Atom in the array
            distance : Distance between atoms i and j
        Returns:
            interaction response

        """
        interaction_type = self.interaction_matrix[j][i]
        interaction_function = self.interaction_map[interaction_type]
        return interaction_function.distance_function[distance]

    def plot_all_interactions(self):
        for interaction in self.interaction_map.values():
            interaction.plot_graph(self.x_axis)
        plt.show()

    def get_initial_coordinates(self):
        atom_coordinates = np.random.uniform(99, 100, (self.number_signals, 3))
        return atom_coordinates

    def get_atom_types(self):
        atom_types = [None for x in self.interaction_matrix]
        while None in atom_types:
            for i in range(len(self.interaction_matrix)):
                for j in range(len(self.interaction_matrix)):
                    if self.interaction_matrix[i][j] == 1:
                        atom_types[i] = "H"
                        atom_types[j] = "H"
                    elif self.interaction_matrix[i][j] == 4:
                        atom_types[i] = "C"
                        atom_types[j] = "C"
                    elif self.interaction_matrix[i][j] == 2:
                        if atom_types[i] == "C":
                            atom_types[j] == "H"
                        elif atom_types[j] == "C":
                            atom_types[i] == "H"
                        elif atom_types[i] == "H":
                            atom_types[j] == "C"
                        elif atom_types[j] == "H":
                            atom_types[i] == "C"
        return atom_types



