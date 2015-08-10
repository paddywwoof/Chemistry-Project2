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
            distance_function[xi] = ((-aad * math.exp(-xi/(10*atc))) + aad + (rad * math.exp(-xi/(10*rtc))) - rad)**power
            if distance_function[xi] > distance_function[xi-1] and distance_function[xi-1] < distance_function[xi-2]:
                print(self.interaction_name+", Minimum at: ", (xi-1)/1000, "Nanometres")
        self.minimum = np.argmin(distance_function)/1000
        self.minimum_y = distance_function[np.argmin(distance_function)]
        self.distance_function = distance_function

    def plot_graph(self, x_axis):
        number_points = 400
        plt.plot(x_axis[:number_points], self.distance_function[:number_points])
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.title('Two Dimensional')
        print("Plotting: " + self.interaction_name)

    def get_response_value(self, distance):
        return self.distance_function[int(distance*1000)]


class DefaultInteraction(Interaction):
    def __init__(self, x_axis, repulsive_amplitude, repulsive_time_constant, depth):
        self.interaction_name = "Default"
        self.depth = depth
        self.distance_function = [0 for i in x_axis]
        for xi in x_axis:
            self.distance_function[xi] = self.depth * (repulsive_amplitude * math.exp(-xi/10*repulsive_time_constant))
        self.minimum = 0.154
        self.minimum_y = self.distance_function[np.argmin(self.distance_function)]


class HMBC(Interaction):
    def __init__(self, interaction_name, x_axis, repulsive_amplitude, repulsive_time_constant, attractive_amplitude,
                 attractive_time_constant, depth, power):
        self.interaction_name = interaction_name
        distance_function = [0 for i in x_axis]

        # To make equation easier to read
        aad = attractive_amplitude * depth
        rad = repulsive_amplitude * depth
        atc = attractive_time_constant
        rtc = repulsive_time_constant
        flat_range = [200, 385]

        xi=0
        while xi < len(x_axis)-1:
            xi += 1
            if xi == flat_range[0]:
                for k in range(*flat_range):
                    distance_function[k] = ((-aad * math.exp(-xi/(10*atc))) + aad + (rad * math.exp(-xi/(10*rtc))) - rad)**power
                xi += (flat_range[1] - flat_range[0])

            distance_function[xi] = ((-aad * math.exp(-xi/(10*atc))) + aad + (rad * math.exp(-xi/(10*rtc))) - rad)**power
            if distance_function[xi] > distance_function[xi-1] and distance_function[xi-1] < distance_function[xi-2]:
                print(self.interaction_name+", Minimum at: ", (xi-1)/1000, "Nanometres")
        self.minimum = np.argmin(distance_function)/1000
        self.minimum_y = distance_function[np.argmin(distance_function)]
        self.distance_function = distance_function
        self.minimum = np.argmin(self.distance_function)/1000
        self.minimum_y = self.distance_function[np.argmin(self.distance_function)]

        print(self.interaction_name+", Minimum at: ", self.minimum, "Nanometres")


class InteractionManager:
    """
    Managers the different types of interaction (e.g. HSQC, COSY etc, and their response values)
    """
    def __init__(self, axis_width, interaction_matrix, type_array, shift_data):
        """
        Args:
            axis_width: Number of points to generate function mapping for

        Returns:
            None
        """
        self.x_axis = range(10*axis_width)
        self.interaction_map = {}
        self.interaction_matrix = interaction_matrix
        self.interaction_matrix_original = np.copy(interaction_matrix)
        self.type_array = type_array
        self.shift_data = shift_data
        self.number_carbon_atoms = self.type_array.count("C")
        self.number_hydrogen_atoms = self.type_array.count("H")
        self.number_atoms = self.number_carbon_atoms + self.number_hydrogen_atoms
        self.bonds, self.bond_orders = self.get_bonds()
        self.bonds = [(x[0]-1, x[1]-1) for x in self.bonds]

        self.global_frag_distances = {
            (9, 10): 0.2467,
            (10, 8): 0.2426,
            (8, 11): 0.2467,
            (11, 9): 0.4945,
            (9, 8): 0.4252,
            (10, 11): 0.4252,
            (7,9): 0.2827,
            (6,9): 0.3041
        }

        for i,j in self.global_frag_distances.keys():
            self.interaction_matrix[i][j] = 6
            self.interaction_matrix[j][i] = 6

    def add_default_interaction(self, index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth):
        """
        Adds default repulsive interaction type
        """
        new_interaction = DefaultInteraction(self.x_axis, repulsive_amplitude, repulsive_time_constant, depth)
        self.interaction_map[index] = new_interaction

    def reset(self):
        self.interaction_matrix = np.copy(self.interaction_matrix_original)

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
               def between(self, x, interval):
        if interval[0] <= x < interval[1]:
            return True
        else:
            return False None
        """
        new_interaction = Interaction(interaction_name, self.x_axis, repulsive_amplitude, repulsive_time_constant, attractive_amplitude,
                                      attractive_time_constant, depth, power)
        self.interaction_map[index] = new_interaction

    def add_hmbc_interaction(self, index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth,
                            attractive_amplitude, attractive_time_constant, power, opt_class = Interaction):
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
               def between(self, x, interval):
        if interval[0] <= x < interval[1]:
            return True
        else:
            return False None
        """
        new_interaction = HMBC(interaction_name, self.x_axis, repulsive_amplitude, repulsive_time_constant, attractive_amplitude,
                                      attractive_time_constant, depth, power)
        self.interaction_map[index] = new_interaction

    def get_bond_length(self, i1, i2):
        """
        Returns Bond Length in Angstroms
        """
        i_type = self.interaction_matrix[i1][i2]
        return self.interaction_map[i_type].minimum

    def calculate_distance(self, v1,v2):
        distance = sqrt( (v1[0]-v2[0]) ** 2 + (v1[1]-v2[1]) ** 2 + (v1[2]-v2[2]) ** 2)
        return min(10000, distance)

    def get_force_coefficient(self, i, j, v1, v2):
        interaction_type = self.interaction_matrix[j][i]

        if interaction_type == 6:
            if (i, j) in self.global_frag_distances:
                scale = self.global_frag_distances[(i, j)]*10
            elif (j, i) in self.global_frag_distances:
                scale = self.global_frag_distances[(j, i)]*10
        else:
            scale = 1
        interaction = self.interaction_map[interaction_type]


        distance = np.linalg.norm(v2-v1)




        delta = (interaction.minimum*scale - distance)

        #print(distance, interaction.minimum*scale, delta)

        return delta*delta



    def interaction_response(self, i, j, v1, v2, debug=False, force=False):
        """
        Args:
            i  : ith Atom index in the array
            j  : jth Atom index in the array
            v1 : Atom i Cartesian Point Vector in nanometres
            v2 : Atom j Cartesian Point Vector in nanometres
        Returns:
            interaction response
        """
        try:
            interaction_type = self.interaction_matrix[j][i]
            if interaction_type == 6:
                if (i, j) in self.global_frag_distances:
                    scale = self.global_frag_distances[(i, j)]*10
                elif (j, i) in self.global_frag_distances:
                    scale = self.global_frag_distances[(j, i)]*10
            else:
                scale = 1
            distance = self.calculate_distance(v1/scale, v2/scale)  # in Nanometers
            interaction = self.interaction_map[interaction_type]
            response_value = interaction.get_response_value(distance)
            if force and False:

                response_value *= -1
                response_value += interaction.minimum_y
            a=abs((response_value - interaction.minimum_y))
            if a > 0.1 and debug and interaction_type!=0:
                print("Interaction %s,%s of type %s is not minimal" % (i+2, j+2, interaction.interaction_name), ": is %s should be %s" % (distance*scale, scale*interaction.minimum), "Response varies by %s "%a)
            elif False:
                print("Interaction %s,%s of type %s is minimal" % (i+2, j+2, interaction.interaction_name), ": is %s should be %s" % (distance*scale, scale*interaction.minimum), "Response varies by %s "%a)
            return response_value
        except IndexError:
            return 0
        except KeyError:
            return 0

    def shape_coordinates(self,atom_coordinates):
        """
        When array is improperly shaped, reshapes the array
        """
        return atom_coordinates.reshape((self.number_atoms, 3))

    def plot_all_interactions(self):

        for interaction in self.interaction_map.values():
            interaction.plot_graph(self.x_axis)
        #pyplot.ion()
        plt.show()

    def get_initial_coordinates(self):
        atom_coordinates = np.random.uniform(99, 101, (self.number_atoms, 3))
        return atom_coordinates

    def between(self, x, interval):
        if interval[0] <= x < interval[1]:
            return True
        else:
            return False

    def get_bonds(self):
        singlebond = [0, 100]
        doublebond = [100, 160]
        oxygen = [160, 222]

        bond_data = set()
        for i in range(len(self.interaction_matrix)):
            for j in range(len(self.interaction_matrix)):
                if self.interaction_matrix[i][j] in [2, 4, 5]:
                    if self.between(self.shift_data[i], singlebond) and self.between(self.shift_data[j], singlebond):
                        bond_order = 1
                    elif self.between(self.shift_data[i], doublebond) and self.between(self.shift_data[j], doublebond):
                        self.interaction_matrix[i][j] = 5
                        bond_order = 2
                    else:
                        bond_order = 1
                    new_bond = [i+1, j+1]
                    new_bond.sort()
                    bond_data.add(tuple(new_bond+[bond_order]))

        bond_data = list(bond_data)
        """
        #Add Oxygens
        for shift_index,shift_value in enumerate(self.shift_data):
            if self.between(shift_value, oxygen):
                self.type_array.append("N")
                bond_data.append([shift_index+1,len(self.type_array),2])
        """
        bonds = [x[0:2] for x in bond_data]
        bond_orders = [x[2] for x in bond_data]

        return bonds, bond_orders

    def print_matrix(self):
        print("     ",*[str(x)+" "*(3-len(str(x))) for x in range(self.number_atoms)])
        for x in range(self.number_atoms):
            print(str(x)+" "*(2-len(str(x))), self.type_array[x], str(self.interaction_matrix[x]).replace(" ", "   "), self.shift_data[x])


