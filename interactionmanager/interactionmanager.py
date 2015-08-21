__author__ = 'martin'
import math
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from filemanager import FileManager
from enum import Enum

class BondLengths:
    single_carbon = 0.154
    double_carbon = 0.142
    carbonyl = 0.123
    single_hydrogen = 0.109

class Interaction:
    def __init__(self, interaction_name, x_axis, repulsive_amplitude, repulsive_time_constant, attractive_amplitude,
                 attractive_time_constant, depth, power):
        self.interaction_name = interaction_name
        distance_function = [0 for i in x_axis]

        # To make equation easier to read
        aa = attractive_amplitude
        ra = repulsive_amplitude
        atc = attractive_time_constant
        rtc = repulsive_time_constant

        for xi in x_axis:
            distance_function[xi] = (depth*(aa * (1 - math.exp(-xi/(10 * atc))) + ra * (math.exp(-xi/(10*rtc))) - 1))**power
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


class BondInteraction(Interaction):
    def __init__(self, interaction_name, bond_length):
        self.interaction_name = interaction_name
        self.minimum = bond_length
        self.minimum_y = 0

    def plot_graph(self, x_axis):
        return

    def get_response_value(self, distance):
        return 0

class DefaultInteraction(Interaction):
    def __init__(self, x_axis, repulsive_amplitude, repulsive_time_constant, depth):
        self.interaction_name = "Default"
        self.depth = depth
        self.distance_function = [0 for i in x_axis]
        for xi in x_axis:
            self.distance_function[xi] = self.depth * (repulsive_amplitude * math.exp(-xi/10*repulsive_time_constant))

            #power=1
            #aa=0

        minimum = 0.154
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

        self.interaction_matrix[self.interaction_matrix == 1] = 0




        self.interaction_matrix_original = np.copy(interaction_matrix)
        self.file_manager = FileManager()
        self.bonds = []
        self.atoms = []
        self.initialise_atoms(type_array, shift_data)

        self.global_frag_distances = {
            (18, 10): 0.3963,
            (18, 15): 0.4834,
            (18, 13): 0.4781,
            (15, 7): 0.3258
        }



        for i, j in self.global_frag_distances.keys():
            self.set_interaction(i, j, 6, False)

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

    def set_interaction(self, x, y, i_type, overwrite=True):
        current_i_type = self.interaction_matrix[x][y]
        if overwrite or current_i_type in [0, 9]:
            self.interaction_matrix[x][y] = i_type
            self.interaction_matrix[y][x] = i_type

    def get_interaction(self, x, y):
        return self.interaction_matrix[x][y]

    def get_all_interaction_atoms(self, interaction_type):
        interactions = []
        for x1, row1 in enumerate(self.interaction_matrix):
            for x2, row2 in enumerate(self.interaction_matrix):
                if x1 < x2:
                    if self.interaction_matrix[x1][x2] == interaction_type:
                        a1 = self.atoms[x1]
                        a2 = self.atoms[x2]
                        interactions.append([a1, a2])  # Gets a list of all HMBC interactions
        return interactions

    def add_bond_interaction(self, index, interaction_name, bond_length):
        new_interaction = BondInteraction(interaction_name, bond_length)
        self.interaction_map[index] = new_interaction

    def add_spatial_interaction(self, index, interaction_name, repulsive_amplitude, repulsive_time_constant, depth,
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

    """
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
        return delta*delta
    """

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
                    scale = self.global_frag_distances[(i, j)]
                elif (j, i) in self.global_frag_distances:
                    scale = self.global_frag_distances[(j, i)]
                else:
                    raise Exception("Scale not defined for %s, %s" % (i, j))
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
                print("Interaction %s,%s of type %s is not minimal" % (i, j, interaction.interaction_name), ": is %s should be %s" % (distance*scale, scale*interaction.minimum), "Response varies by %s "%a)
            elif False:
                print("Interaction %s,%s of type %s is minimal" % (i, j, interaction.interaction_name), ": is %s should be %s" % (distance*scale, scale*interaction.minimum), "Response varies by %s "%a)
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

    def add_atom(self, shift_value, atom_type):
        index_value = len(self.atoms)
        new_atom = Atom(index_value, shift_value, atom_type)
        self.atoms.append(new_atom)

    def initialise_atoms(self, type_array, shift_data):
        for index, atom_type in enumerate(type_array):
            self.add_atom(shift_data[index], atom_type)

        for shift_index, shift_value in enumerate(shift_data):  # Add Inferred C=O Oxygen atoms
            if between(shift_value, oxygen_bond):
                self.add_atom(shift_value, "O")
                self.add_matrix_entry()
                self.set_interaction(shift_index, len(self.atoms)-1, 7)

        bond_set = set()
        for i in range(len(self.interaction_matrix)):
            for j in range(len(self.interaction_matrix)):
                i_type = self.get_interaction(i, j)
                if i_type in [2, 4, 5, 7]:
                    if i_type in [5, 7]:
                        bond_order = 1
                    else:
                        bond_order = 1
                    bond_set.add((self.atoms[min(i, j)], self.atoms[max(i, j)], bond_order, i_type))


        for x in list(bond_set):
            self.add_bond(*x)

    def add_matrix_entry(self):
        im = np.zeros(shape=(len(self.interaction_matrix)+1, len(self.interaction_matrix)+1), dtype=self.interaction_matrix.dtype)
        im[:-1, :-1] = self.interaction_matrix
        self.interaction_matrix = im

    def print_matrix(self):
        print("     ", *[str(x)+" "*(3-len(str(x))) for x in range(self.get_number_atoms())])

        for x, atom in enumerate(self.atoms):
            print(str(x)+" "*(2-len(str(x))), self.atoms[x].atom_type, str(self.interaction_matrix[x]).replace(" ", "   "), self.atoms[x].shift_value)

    def write_numpy_to_mol(self, path, coordinates):
        self.file_manager.write_numpy_to_mol(path, self.bonds, self.atoms, coordinates)

    def add_bond(self, atom1, atom2, bond_order, interaction_type):
        new_bond = Bond(atom1, atom2, bond_order, interaction_type)
        self.bonds.append(new_bond)

    def get_bond_order(self, index):
        return [bond.bond_order for bond in self.bonds]

    def get_bonds_array(self):
        return [bond.get_indices() for bond in self.bonds]

    def get_total_mass(self):
        total_mass = 0
        if len(self.atoms) == 0:
            raise Exception("Atom list is empty")
        for atom in self.atoms:
            total_mass += mass_dict[atom.atom_type]
        return total_mass

    def get_type(self, index):
        return self.atoms[index].atom_type

    def get_type_array(self):
        return [atom.atom_type for atom in self.atoms]

    def get_free_valencies(self):
        return [atom.get_free_valency() for atom in self.atoms]


    """
    def update_bond(self, atom1, atom2, bond_order):
        atom1.update_bond(atom2.index_value, bond_order)
        atom2.update_bond(atom1.index_value, bond_order)
        for i, bond in enumerate(self.bonds):
            if bond[0:2] in [[atom1.index_value, atom2.index_value], [atom1.index_value, atom2.index_value]]:
                self.bonds[i][2] = bond_order
    """

    def update_all_bonds(self):
        # Keeps restarting updating bonds until bonds don't update any more
        for bond in self.bonds:
            if bond.update_bond():
                self.update_all_bonds()
                return

    """
    def check_double_bond_rings(self, atom_chain):
        #  Checks for rings with all double bond shifts (usually benzene)
        end_chain_atom = atom_chain[-1]
        for atom in [self.atoms[x] for x in end_chain_atom.get_adjacent()]:
            if atom.needs_double_bond():
                if atom not in atom_chain:
                    self.check_double_bond_rings(atom_chain+[atom])
                elif atom != atom_chain[-2]:
                    self.update_bond(atom[-1], atom[-2], 2)
    """

atom_valencies = {"C": 4, "H": 1, "O": 2}
singlebond = [0, 100]
doublebond = [100, 160]
oxygen_bond = [160, 222]
mass_dict = {"C": 12, "O": 16, "H": 1, "N": 14}

def between(x, interval):
    if interval[0] <= x < interval[1]:
        return True
    else:
        return False


class Atom:
    def __init__(self, index_value, shift_value, atom_type):
        self.index_value = index_value
        self.shift_value = shift_value
        self.atom_type = atom_type
        self.total_valency = atom_valencies[atom_type]
        self.bonds = []

        if between(self.shift_value, doublebond):
            self.has_double_bond = True
        else:
            self.has_double_bond = False

    def needs_double_bond(self):
        if self.has_double_bond and 2 not in [x.bond_order for x in self.bonds]:
            if self.atom_type == "H":
                raise Exception("Hydrogen Cannot Double Bond")
            return True
        else:
            return False

    def get_free_valency(self):
        valency = self.total_valency
        for bond in self.bonds:
            valency -= bond.bond_order
        if self.needs_double_bond():
            valency -= 1
        if valency < 0:
            raise Exception("Negative Valency Error: %s Atom %s has negative valency" %(self.atom_type, self.index_value))
        return valency

    def __str__(self):
        return str(self.index_value)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return self.index_value

    def get_adjacent(self):
        return [bond.other(self) for bond in self.bonds]

    def add_bond(self, bond):
        if bond not in self.bonds:
            self.bonds.append(bond)

    def get_bond(self, atom):
        for bond in self.atoms:
            if atom in bond:
                return bond

class Bond:
    def __init__(self, atom1, atom2, bond_order, inferred_by, aromatic=False):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_indices = [atom1.index_value, atom2.index_value]
        self.bond_order = 1
        self.bond_length = 10
        self.inferred_by = inferred_by
        self.aromatic = aromatic
        if bond_order != 2 and aromatic:
            raise Exception("Only double bonds can be aromatic")
        self.atom1.add_bond(self)
        self.atom2.add_bond(self)

    def __getitem__(self, key):
        return self.bond_indices[key]

    def update_bond(self):
        if "H" in [self.atom1.atom_type, self.atom2.atom_type]:
            if self.bond_length != BondLengths.single_hydrogen:
                self.bond_length = BondLengths.single_hydrogen
                return True

        elif "O" in [self.atom1.atom_type, self.atom2.atom_type]:  # If C<>O Bond
            if self.atom1.atom_type == "O":
                oxygen = self.atom1
            else:
                oxygen = self.atom2
            if between(oxygen.shift_value, oxygen_bond):
                if self.bond_order != 2:
                    self.bond_order = 2
                    self.bond_length = BondLengths.carbonyl
                    return True

        elif ["C", "C"] == [self.atom1.atom_type, self.atom2.atom_type]:
            if self.atom1.needs_double_bond() and self.atom2.needs_double_bond():
                """
                Atoms bonded to atom1 that need double bond and vice verse
                If atom1 has only one double bond candidate and it is atom2, or vice versa
                then atom1-atom2 is a double bond
                """
                atom1_doublebond_adjacents = [atom for atom in self.atom1.get_adjacent() if atom.needs_double_bond()]
                atom2_doublebond_adjacents = [atom for atom in self.atom2.get_adjacent() if atom.needs_double_bond()]

                a = len(atom1_doublebond_adjacents) == 1 and atom1_doublebond_adjacents[0] == self.atom2
                b = len(atom2_doublebond_adjacents) == 1 and atom2_doublebond_adjacents[0] == self.atom1
                if (a or b) and (self.atom1.get_free_valency() == 0 or self.atom2.get_free_valency() == 0):
                    if "H" in [self.atom1.atom_type, self.atom2.atom_type]:
                        raise Exception("Hydrogen Cannot Double Bond")
                    if self.bond_order != 2:
                        self.bond_order = 2
                        self.bond_length = BondLengths.double_carbon
                        return True

            else:
                if self.bond_length != BondLengths.single_carbon:
                    self.bond_length = BondLengths.single_carbon
                    return True

        return False

    def check_double_bond_rings(self, atom_chain=[]):
        #  Checks for rings with all double bond shifts (usually benzene)
        if len(atom_chain) == 0:
            atom_chain = [self.atom1]
        end_chain_atom = atom_chain[-1]
        for atom in end_chain_atom.get_adjacent():
            if atom.needs_double_bond():
                if atom not in atom_chain:
                    return self.check_double_bond_rings(atom_chain+[atom])
                elif atom != atom_chain[-2]:
                    for x in atom_chain:
                        if len([y for y in x.get_adjacent() if y in atom_chain]) > 2:
                            return
                    atom_chain = atom_chain[:-1]
                    if len(atom_chain) % 2 != 0:
                        raise Exception("This Double Bond Ring Should not be possible")
                    while len(atom_chain) > 0:
                        bond = atom_chain[-1].get_bond(atom_chain[-2])
                        bond.bond_order = 2
                        bond.aromatic = True
                        atom_chain = atom_chain[:-2]
                    return True
        return False


    def other(self, atom):
        if atom == self.atom1:
            return self.atom2
        elif atom == self.atom2:
            return self.atom1
        else:
            raise Exception("Atom not in bond")

    def __str__(self):
        return str(self.get_indices())

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return int("".join([str(ord(x)) for x in self.__str__()]))

    def __contains__(self, item):
        return item in [self.atom1, self.atom2]

    def get_indices(self):
        a1 = min([self.atom1.index_value, self.atom2.index_value])
        a2 = max([self.atom1.index_value, self.atom2.index_value])
        return [a1, a2]

