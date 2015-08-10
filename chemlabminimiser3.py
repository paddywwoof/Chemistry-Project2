from filemanager import FileManager
import numpy as np
import random
from chemlab.graphics.transformations import rotation_matrix
from graphics import MolecularGraphics

from interactionmanager import InteractionManager
from signalmanager import OneDSignalManager, TwoDSignalManager

class ChemLabMinimiser:
    def __init__(self):
        self.interaction_manager = get_interaction_manager(*get_twod_signal_manager().get_interaction_data())
        # Consider 4-Bond Interactions
        im = self.interaction_manager.interaction_matrix
        self.interaction_matrix = self.interaction_manager.interaction_matrix
        self.type_array = list(self.interaction_manager.type_array)
        self.shift_data = list(self.interaction_manager.shift_data)
        self.number_atoms = self.interaction_manager.number_atoms
        self.bonds = list(self.interaction_manager.bonds)
        self.bond_orders = list(self.interaction_manager.bond_orders)
        self.best_response_value = 1000000.0
        ac = np.random.rand(self.number_atoms, 3)
        #1ac = FileManager().read_numpy_from_xyz('manual2.xyz') * 0.1
        self.coordinate_manager = CoordinateManager(self, ac, self.interaction_manager)
        self.iteration_number = 0
        self.fragments = self.generate_fragments_list()
        self.graphics = MolecularGraphics(self.coordinate_manager.get_coordinates(), self.type_array, np.array(self.bonds))

    def calculate_valencies(self):
        valencies = {"H":1, "C":4, "N":4}
        self.free_valencies = [valencies[x] for x in self.interaction_manager.type_array]
        for bond in self.bonds:
            self.free_valencies[bond[0]]-=1
            self.free_valencies[bond[1]]-=1
        for i, shift in enumerate(self.shift_data):
            if 100 < shift < 160:
                self.free_valencies[i]-=1
        """
        for index, type in enumerate(self.type_array):
            if self.free_valencies[index] == 0 and self.type_array[index] == "C":
                self.type_array[index]="N"
        """
        #input(self.free_valencies)


    #TODO: Merge generate functions together?
    def generate_fragments_list(self):
        """
        Creates a list of fragment objects
        """
        self.calculate_valencies()
        fragments = []
        fragment_data = self.generate_fragment_groups()
        for fragment_piece in fragment_data:
            global_indices = fragment_piece['global indices']
            global_bond_indices = fragment_piece['global bond indices']
            bond_orders = fragment_piece['bond orders']
            fragment = Fragment(global_indices,
                                global_bond_indices,
                                self.coordinate_manager,
                                self.interaction_manager,
                                self.free_valencies
                                )
            fragments.append(fragment)
        return fragments

    def generate_fragment_groups(self):
        fragment_groups = []
        unselected = [x for x in range(self.number_atoms)]
        while len(unselected) > 0:
            new_atoms = set()
            new_atoms.add(unselected[0])
            atoms = []
            new_bonds = []
            new_bond_orders = []
            while len(new_atoms) > 0:
                a = list(new_atoms)[0]
                atoms.append(a)
                new_atoms.remove(a)
                for bond_index, bond in enumerate(self.bonds):
                    if atoms[-1] == bond[0] and bond[1] not in atoms:
                        new_atoms.add(bond[1])
                    if atoms[-1] == bond[1] and bond[0] not in atoms:
                        new_atoms.add(bond[0])
                    if a in bond and list(bond) not in new_bonds:
                        new_bonds.append(list(bond))
                        new_bond_orders.append(self.bond_orders[bond_index])
            atoms.sort()
            fragment_piece = dict()
            fragment_piece["global indices"] = atoms
            fragment_piece["global bond indices"] = new_bonds
            fragment_piece["bond orders"] = new_bond_orders
            fragment_groups.append(fragment_piece)
            unselected = [x for x in range(self.number_atoms) if x not in sum([y["global indices"] for y in fragment_groups],[])]
        return fragment_groups

    def main(self):
        print(self.coordinate_manager.calculate_response_value(True))
        self.interaction_manager.print_matrix()
        input("")

        print(self.coordinate_manager.calculate_response_value())

        for x in range(0,10):
            self.iteration2()

        self.graphics.run(self.iteration)

    def iteration2(self):
        self.calculate_valencies()
        atom_coordinates = self.coordinate_manager.get_coordinates()
        hmbcs = []
        for x1, row1 in enumerate(self.interaction_matrix):
            for x2, row2 in enumerate(self.interaction_matrix):
                if x1 < x2:
                    if self.interaction_matrix[x1][x2] == 3:
                        hmbcs.append([x1, x2]) # Gets a list of all HMBC interactions

        hmbcs = [x for x in hmbcs if self.bonddistance(x[0], x[1]) not in [2, 3]]


        for hmbc in hmbcs: #Iterate though each interaction
            for frag1 in self.fragments:
                for frag2 in self.fragments:
                    if frag1 != frag2 or len(self.fragments) == 1:
                        if not(hmbc[0] in frag1.global_indices and hmbc[1] in frag2.global_indices) and not(hmbc[1] in frag1.global_indices and hmbc[0] in frag2.global_indices):
                            continue
                        if self.free_valencies[max(hmbc)]==0:
                            continue

                        if self.type_array[min(hmbc)] == "H":
                            bonded_carbon_bonds = [x for x in self.bonds if min(hmbc) in x] # Index of carbon bonded to hydrogen
                            if len(bonded_carbon_bonds)>1:
                                raise Exception("Something Wrong Here")
                            else:
                                bonded_carbon = max(bonded_carbon_bonds[0])



                            if self.free_valencies[bonded_carbon] != 0: # If Bonded carbon can be bonded to, ignore
                                #print("Carbon %s can be bonded to still!"%(bonded_carbon+2))
                                continue
                            else:
                                #print("Carbon %s cannot be bonded to, find neighbours!"%(bonded_carbon+2))

                                # Look for carbons adjacent to this carbon and check if they have free valencies
                                # If there is only one, then this must be the correct bond

                                adjacent_carbon_bonds = [list(x) for x in self.bonds if bonded_carbon in x]

                                three_bond_carbons = [x for x in sum(adjacent_carbon_bonds,[]) if x!=bonded_carbon]
                                adjacent_carbons = [x for x in list(set(three_bond_carbons)) if self.free_valencies[x] != 0 and self.type_array[x] in ["N","C"]]

                                if len(adjacent_carbons) == 1:
                                    bond = [adjacent_carbons[0], max(hmbc)]
                                    if bond[0]==bond[1] or bond in self.bonds:
                                        continue
                                    print(bonded_carbon)
                                    print("HMBC:", hmbc)
                                    print("Bond Chain: ", min(hmbc), bonded_carbon, adjacent_carbons[0], max(hmbc))
                                    print(adjacent_carbons)
                                    self.addbond(frag1, frag2, bond)
                                    return

    def bonddistance(self,i,j):
        fragment = [x for x in self.fragments if i in x.global_indices][0]

        one_bond_atoms = self.get_adjacent(fragment, i)
        two_bond_atoms = self.uniq([self.get_adjacent(fragment, x, one_bond_atoms+[i]) for x in one_bond_atoms])
        three_bond_atoms = self.uniq([self.get_adjacent(fragment, x, one_bond_atoms+two_bond_atoms+[i]) for x in two_bond_atoms])
        if i == j:
            return 0
        if j in one_bond_atoms:
            return 1
        if j in two_bond_atoms:
            return 2
        if j in three_bond_atoms:
            return 3

    def get_adjacent(self, fragment, i, expanded=[]):
        adjacent_bonds = [x for x in fragment.global_bond_indices if i in x]
        adjacent_atoms = [x for x in list(set(sum(adjacent_bonds, []))) if x != i and x not in expanded]
        return adjacent_atoms

    def uniq(self,x):
        return list(set(sum(x,[])))















    def iteration(self):
        print("Starting Iteration %s with %s fragments"%(self.iteration_number,len(self.fragments)))
        self.iteration_number += 1
        self.graphics.update_coordinates(self.coordinate_manager.get_coordinates())
        if self.iteration_number < 1000:
            self.optimise_substructures()
        else:
            self.form_new_bonds()
            self.iteration_number = 0


    def optimise_substructures(self):
        for fragment in self.fragments:
            fragment.start_iter()
            fragment.rotate_bonds(sequence=1)
            fragment.rotate_bonds(sequence=1)
            fragment.verify()
            fragment.start_iter()
            fragment.translate()
            fragment.verify()
            self.graphics.update_coordinates(self.coordinate_manager.get_coordinates())

    def form_new_bonds(self):
        shortest_bond = None
        closest_fragment = None
        shortest_distance = 1000
        atom_coordinates = self.coordinate_manager.get_coordinates()

        fragment1 = random.choice(self.fragments)
        for fragment2 in self.fragments:
            if fragment1 != fragment2 or len(self.fragments)==1:
                for atom1_index in fragment1.global_indices:
                    for atom2_index in fragment2.global_indices:
                        frag1_atom = atom_coordinates[atom1_index]
                        frag2_atom = atom_coordinates[atom2_index]
                        if np.linalg.norm(frag2_atom-frag1_atom) < shortest_distance and atom2_index!=atom1_index:
                            if self.free_valencies[atom1_index] > 0 and self.free_valencies[atom2_index] > 0:
                                a = [not(fragment1.has_hydrogens[atom1_index] and fragment1.has_hydrogens[atom2_index]), not([atom1_index, atom2_index] in self.bonds), not([atom2_index, atom1_index] in self.bonds)]
                                if False not in a:
                                    shortest_bond = [atom1_index, atom2_index]
                                    shortest_distance = np.linalg.norm(frag2_atom-frag1_atom)
                                    closest_fragment = fragment2
                                else:
                                    print(self.bonds)
                                    print(atom1_index, atom2_index)
                                    print(a)
        print("For Fragment with atoms %s"%fragment1.global_indices, "New Bond %s"%str(shortest_bond))
        if closest_fragment == None:
            return
        else:
            self.coordinate_manager.reset()
        self.addbond(fragment2, closest_fragment, shortest_bond)

    def addbond(self,fragment1, fragment2,bond):

        self.interaction_manager.interaction_matrix[bond[0]][bond[1]] = 4
        self.interaction_manager.interaction_matrix[bond[1]][bond[0]] = 4

        print("Adding Bond %s"%str(bond))
        self.bonds.append(list(bond))
        self.calculate_valencies()

        if len(self.fragments)>1:
            self.merge(fragment1, fragment2, bond)
            self.fragments.remove(fragment1)
            if fragment1!=fragment2:
                self.fragments.remove(fragment2)


        self.coordinate_manager.reset()
        self.graphics.update(self.coordinate_manager.get_coordinates(), self.type_array, self.bonds)


    def merge(self, fragment1, fragment2, shortest_bond):
        """
        Creates a list of fragment objects
        """
        new_fragment_indices = fragment1.global_indices + fragment2.global_indices
        new_fragment_bonds = fragment1.global_bond_indices + fragment2.global_bond_indices + [shortest_bond]
        merged_fragment = Fragment(new_fragment_indices, new_fragment_bonds, self.coordinate_manager, self.interaction_manager, self.free_valencies)
        self.fragments.append(merged_fragment)
        return merged_fragment

    def checkHMBC(self, i1, i2, frag1, frag2):
        return True

        all_bonds = frag1.global_bond_indices + frag2.global_bond_indices + [[i1,i2]]

        all_indices = list(set(sum(all_bonds,[])))

        i1_hmbcs = [x for x in range(0,self.number_atoms) if self.interaction_matrix[i1][x]==3 and x in all_indices]

        one_bond_gap = []
        for bond in all_bonds:
            if i1 in bond:
                if i1 == bond[0]:
                    one_bond_gap.append(bond[1])
                else:
                    one_bond_gap.append(bond[0])

        two_bond_gap = []
        for index in one_bond_gap:
            for bond in all_bonds:
                if index in bond and index not in one_bond_gap and index != i1:
                    if index == bond[0]:
                        two_bond_gap.append(bond[1])
                    else:
                        two_bond_gap.append(bond[0])

        three_bond_gap = []
        for index in two_bond_gap:
            for bond in all_bonds:
                if index in bond and index not in one_bond_gap and index not in two_bond_gap and index != i1:
                    if index == bond[0]:
                        three_bond_gap.append(bond[1])
                    else:
                        three_bond_gap.append(bond[0])

        delta_indices =  [x for x in i1_hmbcs if x not in two_bond_gap+three_bond_gap]
        if len(delta_indices) > 0:
            print(delta_indices)
            return False
        else:
            return True
        i2_hmbcs = [x for x in range(0,self.number_atoms) if self.interaction_matrix[i2][x]==3 and x in all_indices]


class Fragment:
    def __init__(self, global_indices, global_bond_indices, coordinate_manager, interaction_manager, free_valencies):
        self.global_indices = global_indices
        self.global_bond_indices = global_bond_indices
        self.type_array = np.array(interaction_manager.type_array)
        self.coordinate_manager = coordinate_manager
        self.interaction_manager = interaction_manager
        #self.project_bond_lengths()
        self.number_atoms = len(self.global_indices)
        self.best_response_value = 1000
        self.verify_index = 0
        self.free_valencies = free_valencies

        self.has_hydrogens = [False for x in range(len(self.coordinate_manager.atom_coordinates))]
        for x in self.global_bond_indices:
            if self.type_array[x[0]]=="H":
                self.has_hydrogens[x[1]]=True
            if self.type_array[x[1]]=="H":
                self.has_hydrogens[x[0]]=True

    def project_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.global_bond_indices:
            if self.interaction_manager.interaction_matrix[bond[0]][bond[1]]==9:
                continue

            frozen, unfrozen = self.bisect_on_bond(bond)
            bond_length = self.interaction_manager.get_bond_length(*bond)
            bond_vector = atom_coordinates[bond[1]] - atom_coordinates[bond[0]]
            corrected_bond_vector = bond_vector * (bond_length/np.linalg.norm(bond_vector))
            atom_coordinates[unfrozen] += (corrected_bond_vector-bond_vector)
        self.coordinate_manager.update(atom_coordinates)

    def reduce_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.global_bond_indices:
            if self.interaction_manager.interaction_matrix[bond[0]][bond[1]]==9:
                continue

            frozen, unfrozen = self.bisect_on_bond(bond)
            bond_length = self.interaction_manager.get_bond_length(*bond)
            bond_vector = atom_coordinates[bond[1]] - atom_coordinates[bond[0]]
            corrected_bond_vector = bond_vector * (bond_length/np.linalg.norm(bond_vector))

            x = random.uniform(0,0.5)
            atom_coordinates[unfrozen] += (corrected_bond_vector-bond_vector) * x
            atom_coordinates[frozen] -= (corrected_bond_vector-bond_vector) * x
        self.coordinate_manager.update(atom_coordinates)


    def bisect_on_bond(self, bond, freeze_left=True):
        """

        """
        bonds = self.global_bond_indices
        frozen_atoms = set()
        new_atoms = set()
        if freeze_left:
            new_atoms.add(bond[0])
            block = bond[1]
        else:
            new_atoms.add(bond[1])
            block = bond[0]
        while len(new_atoms) > 0:
            a = list(new_atoms)[0]

            for bond in bonds:
                if block not in bond:
                    if a == bond[0] and a not in frozen_atoms:
                        new_atoms.add(bond[1])
                    if a == bond[1] and a not in frozen_atoms:
                        new_atoms.add(bond[0])
            frozen_atoms.add(a)
            new_atoms.remove(a)
        unfrozen_atoms = [x for x in list(set(sum(bonds, []))) if x not in frozen_atoms]#


        return list(frozen_atoms), unfrozen_atoms

    def rotate_bonds(self, bond_index=None, angle=None, sequence=1):
        self.reduce_bond_lengths()
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        if len(self.global_bond_indices) > 0:
            bond = random.choice(self.global_bond_indices)
            if random.randrange(0,2)==0:
                bond.reverse()
            index1 = bond[0]
            index2 = bond[1]
            point_a = np.copy(new_atom_coordinates[index1])
            point_b = np.copy(new_atom_coordinates[index2])
            uf = self.bisect_on_bond(bond, True)[1]
            angle = (np.pi/100) * random.uniform(0, 25)

            if random.randrange(0,2)==0:
                rvec = np.random.rand(3)
                rotation_axis = np.cross(point_b-point_a, rvec)
            else:
                rotation_axis = point_b - point_a

            rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
            new_atom_coordinates -= point_b
            new_atom_coordinates[uf] = np.dot(new_atom_coordinates[uf], rotation.T)
            new_atom_coordinates += point_b
            self.coordinate_manager.update(new_atom_coordinates)
        else:
            self.translate()

    def verify(self):
        return self.coordinate_manager.verify_global()

    def start_iter(self):
        self.coordinate_manager.start_iter()

    def translate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()

        xvec = np.random.normal(loc=1.0,scale=3.0)*(np.random.rand(3)-0.5)
        xvec /= np.linalg.norm(xvec)

        vec = random.uniform(0,0.6)*xvec
        new_atom_coordinates[self.global_indices] += vec

        self.coordinate_manager.update(new_atom_coordinates)

    def rotate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        atom_index = random.choice(self.global_indices)
        atom = new_atom_coordinates[atom_index]
        new_atom_coordinates -= atom
        angle = np.random.normal(np.pi*0.1, np.pi*0.2 )
        rotation_axis = np.random.rand(3)
        rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
        new_atom_coordinates[self.global_indices] = np.dot(new_atom_coordinates[self.global_indices], rotation.T)
        new_atom_coordinates += atom
        self.coordinate_manager.update(new_atom_coordinates)



class CoordinateManager:
    def __init__(self, parent, atom_coordinates, interaction_manager):
        self.parent = parent
        self.atom_coordinates = atom_coordinates
        self.old_atom_coordinates = np.copy(atom_coordinates)
        self.interaction_manager = interaction_manager
        self.best_response_value = self.calculate_response_value()
        self.file_manager = FileManager()
        self.temperature = 10
        self.number_atoms = self.interaction_manager.number_atoms

    def start_iter(self):
        self.old_atom_coordinates = np.copy(self.atom_coordinates)

    def reset(self):
        self.best_response_value = 1000000

    def update(self, atom_coordinates):
        for x in [0, 1, 2]:
            atom_coordinates[:, x] -= min(atom_coordinates[:, x])
        self.atom_coordinates = atom_coordinates

    def verify_global(self):
        response = self.calculate_response_value()
        if response < self.best_response_value:
            self.best_response_value = response
            if False:
                print("New Best Response: %s at iteration %s"%(response,self.parent.iteration_number))
            self.file_manager.write_numpy_to_mol("tempfile.mol", self.parent.bonds, self.parent.type_array, self.atom_coordinates)
            return True
        else:
            self.revert()
            return False

    def revert(self):
        self.atom_coordinates = self.old_atom_coordinates

    def calculate_response_value(self, debug=False):
        response = 0
        for i, v1 in enumerate(self.atom_coordinates):
            for j, v2 in enumerate(self.atom_coordinates):
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, v1, v2, debug)
        return response

    def get_coordinates(self):
        return np.copy(self.atom_coordinates)


def get_twod_signal_manager():
    oned_signal_manager = OneDSignalManager()
    oned_signal_manager.add_nmr_signals('resources/nmr/oned/hydrogen_integration_data.txt', "H")
    oned_signal_manager.add_nmr_signals('resources/nmr/oned/carbon_integration_data.txt', "C")
    twod_signal_manager = TwoDSignalManager(oned_signal_manager)
    twod_signal_manager.add_nmr_signals("COSY", 'resources/nmr/twod/cosy/cosy_peak_data.txt')
    twod_signal_manager.add_nmr_signals("HSQC", 'resources/nmr/twod/hsqc/hsqc_peak_data.txt')
    twod_signal_manager.add_nmr_signals("HMBC", 'resources/nmr/twod/hmbc/hmbc_peak_data.txt')
    return twod_signal_manager

def get_interaction_manager(interaction_matrix, type_array, shift_data):
    interaction_manager = InteractionManager(1001, interaction_matrix, type_array, shift_data)
    interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=0.2, depth=10)
    interaction_manager.add_new_interaction(index=1, interaction_name="COSY 3-Bond H-H          ", repulsive_amplitude=0.8, repulsive_time_constant=6.28, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=2, interaction_name="HSQC 1-Bond H-C          ", repulsive_amplitude=0.8, repulsive_time_constant=2.26, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_hmbc_interaction(index=3, interaction_name="HMBC 2/3 Bond H-C        ", repulsive_amplitude=0.8, repulsive_time_constant=6.50, depth=3.8, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=4, interaction_name="INAD Single 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.49, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=5, interaction_name="INAD Double 1-Bond C=C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.15, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=6, interaction_name="NOESY Normalised Function", repulsive_amplitude=0.8, repulsive_time_constant=2.03, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    return interaction_manager



def end():
    input("Press Any Key to Quit")
    raise SystemExit

def main():
    clm = ChemLabMinimiser()
    file_manager = FileManager()
    clm.main()
    end()

if __name__ == "__main__":
    #import cProfile
    #cProfile.run('main()')
    main()
