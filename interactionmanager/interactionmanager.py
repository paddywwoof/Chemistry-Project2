__author__ = 'martin'
import math
from math import sqrt

import numpy as np
import matplotlib.pyplot as plt
from filemanager import FileManager


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
            distance_function[xi] = (aad * (1 - math.exp(-xi/(10 * atc))) + rad * (math.exp(-xi/(10*rtc))) - 1)**power
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
        self.file_manager = FileManager()



        self.type_array = type_array
        self.shift_data = shift_data
        self.bonds, self.bond_orders, bond_data = self.generate_bonds()
        self.bonds = [(x[0], x[1]) for x in self.bonds]

        self.atoms = []
        for i in range(len(self.type_array)):
            new_atom = Atom(shift_data[i],  type_array[i], [bond for bond in bond_data if i in bond[0:2]])
            self.atoms.append(new_atom)




        self.global_frag_distances = {
        (2 ,3 ): 1.767,
        (2 ,4 ): 1.769,
        (2 ,15): 2.431,
        (2 ,14): 2.500,
        (2 ,17): 4.878,
        (2 ,13): 4.133,
        (3 ,15): 3.060,
        (3 ,14): 2.486,
        (3 ,13): 2.638,
        (3 ,12): 4.110,
        (3 ,11): 3.835,
        (4 ,14): 3.072,
        (4 ,15): 2.526,
        (4 ,16): 4.788,
        (4 ,17): 4.464,
        (4 ,11): 4.522,
        (4 ,12): 4.084,
        (4 ,13): 2.896,
        (14,15): 1.763,
        (14,13): 4.101,
        (14,11): 4.799,
        (14,12): 4.722,
        (14,20): 3.555,
        (14,16): 4.505,
        (14,17): 3.594,
        (14,6 ): 4.864,
        (15,13): 4.359,
        (15,12): 4.864,
        (15,17): 2.473,
        (15,16): 3.583,
        (15,6 ): 4.523,
        (15,7 ): 4.840,
        (16,17): 1.763,
        (16,12): 4.162,
        (16,20): 4.649,
        (16,5 ): 2.441,
        (16,6 ): 3.069,
        (16,7 ): 2.531,
        (17,5 ): 3.046,
        (17,6 ): 2.562,
        (17,7 ): 2.391,
        (5 , 6): 1.780,
        (5 , 7): 1.759,
        (5 ,18): 3.128,
        (5 ,19): 4.337,
        (5 ,10): 4.550,
        (5 , 9): 3.398,
        (5 , 8): 4.741,
        (5 ,20): 4.863,
        (6 ,20): 4.783,
        (6 ,19): 3.615,
        (6 , 9): 4.801,
        (7 ,18): 4.566,
        (7 , 9): 4.866,
        (20,13): 4.076,
        (20,11): 3.334,
        (20,12): 3.675,
        (20,19): 4.772,
        (13,11): 1.769,
        (13,12): 1.776,
        (12,11): 1.772,
        (18,19): 1.779,
        (19,10): 2.486,
        (19,8 ): 2.498,
        (19,9 ): 3.068,
        (18,9 ): 2.509,
        (18,10): 3.064,
        (18,8 ): 2.483
        }

        for i,j in self.global_frag_distances.keys():
            self.interaction_matrix[i-2][j-2] = 6
            self.interaction_matrix[j-2][i-2] = 6

    def add_default_interaction(self, index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth):
        """
        Adds default repulsive interaction type
        """
        new_interaction = DefaultInteraction(self.x_axis, repulsive_amplitude, repulsive_time_constant, depth)
        self.interaction_map[index] = new_interaction

    def get_number_atoms(self):
        return len(self.interaction_matrix)

    def reset(self):
        self.interaction_matrix = np.copy(self.interaction_matrix_original)

    def set_interaction(self, x, y, i_type):
        self.interaction_matrix[x][y] = i_type
        self.interaction_matrix[y][x] = i_type

    def get_interaction(self,x,y):
        return self.interaction_matrix[x][y]


    def get_all_interactions(self, interaction_type):
        interactions = []
        for x1, row1 in enumerate(self.interaction_matrix):
            for x2, row2 in enumerate(self.interaction_matrix):
                if x1 < x2:
                    if self.interaction_matrix[x1][x2] == interaction_type:
                        interactions.append([x1, x2])  # Gets a list of all HMBC interactions
        return interactions

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

    def calculate_distance(self, v1, v2):
        distance = sqrt( (v1[0]-v2[0]) ** 2 + (v1[1]-v2[1]) ** 2 + (v1[2]-v2[2]) ** 2)
        return min(10000, distance)

    def get_force_coefficient(self, i, j, v1, v2):
        interaction_type = self.interaction_matrix[j][i]

        if interaction_type == 6:
            if (i+2, j+2) in self.global_frag_distances:
                scale = self.global_frag_distances[(i+2, j+2)]
            elif (j+2, i+2) in self.global_frag_distances:
                scale = self.global_frag_distances[(j+2, i+2)]
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
                if (i+2, j+2) in self.global_frag_distances:
                    scale = self.global_frag_distances[(i+2, j+2)]
                elif (j+2, i+2) in self.global_frag_distances:
                    scale = self.global_frag_distances[(j+2, i+2)]
                else:
                    raise Exception("Scale not defined for %s, %s"%(i,j))
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

    def plot_all_interactions(self):

        for interaction in self.interaction_map.values():
            interaction.plot_graph(self.x_axis)
        #pyplot.ion()
        plt.show()

    def between(self, x, interval):
        if interval[0] <= x < interval[1]:
            return True
        else:
            return False

    def generate_bonds(self):
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
                        self.set_interaction(i, j, 5)
                        bond_order = 2
                    else:
                        bond_order = 1
                    new_bond = [i, j]
                    new_bond.sort()
                    bond_data.add(tuple(new_bond+[bond_order]))

        bond_data = list(bond_data)

        #Add Oxygens
        for shift_index,shift_value in enumerate(self.shift_data):
            if self.between(shift_value, oxygen):
                self.type_array.append("O")
                self.shift_data.append(0)
                bond_data.append([shift_index,len(self.type_array)-1 ,2])
                im = np.zeros(shape=(len(self.interaction_matrix)+1, len(self.interaction_matrix)+1), dtype=self.interaction_matrix.dtype)
                im[:-1, :-1] = self.interaction_matrix
                self.interaction_matrix = im
                self.set_interaction(shift_index, len(self.type_array)-1, 5)

        bonds = [x[0:2] for x in bond_data]

        bond_orders = [x[2] for x in bond_data]
        return bonds, bond_orders, bond_data

    def print_matrix(self):
        print("     ",*[str(x)+" "*(3-len(str(x))) for x in range(self.get_number_atoms())])
        for x in range(self.get_number_atoms()):
            self.type_array[x]
            self.interaction_matrix[x]
            self.shift_data[x]
            print(str(x)+" "*(2-len(str(x))), self.type_array[x], str(self.interaction_matrix[x]).replace(" ", "   "), self.shift_data[x])

    def get_total_mass(self):
        total_mass = 0
        mass_dict = {"C":12, "O":16, "H":1, "N":14}
        for x in self.type_array:
            total_mass += mass_dict[x]
        return total_mass

    def add_bond(self, bond, bond_order):
        self.bonds.append(tuple(bond))
        self.bond_orders.append(bond_order)

    def get_type(self, index):
        return self.type_array[index]

    def get_bond_order(self, index):
        return self.bond_orders[index]

    def get_bonds(self):
        return self.bonds

    def get_type_array(self):
        return self.type_array

    def write_numpy_to_mol(self, path, coordinates):
        self.file_manager.write_numpy_to_mol(path, self.bonds,self.bond_orders, self.type_array, coordinates)

    def get_free_valencies(self):
        valencies = {"H":1, "C":4, "N":4, "O":2}
        free_valencies = [valencies[x] for x in self.type_array]
        for bond in self.bonds:
            free_valencies[bond[0]] -= 1
            free_valencies[bond[1]] -= 1
        for i, shift in enumerate(self.shift_data):
            if 100 < shift < 160:
                free_valencies[i] -= 1
        return free_valencies



atom_valencies = {"C":4, "H":1}
singlebond = [0, 100]
doublebond = [100, 160]
oxygen = [160, 222]

def between(x, interval):
    if interval[0] <= x < interval[1]:
        return True
    else:
        return False


class Atom:
    def __init__(self, shift_value, atom_type, bonds):
        self.shift_value = shift_value
        self.atom_type = atom_type
        self.total_valency = atom_valencies[atom_type]
        self.bonds = bonds

        if between(self.shift_value, doublebond):
            self.has_double_bond = True
        else:
            self.has_double_bond = False

    def get_free_valency(self):
        valency = self.total_valency
        for bond in self.bonds:
            valency -= bond[2]
        if self.has_double_bond and 2 not in [x[2] for x in self.bonds]:
            valency -= 1

    def get_adjacent(self):
        pass

    def add_bond(self):
        pass


