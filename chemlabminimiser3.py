import numpy as np
import random
from chemlab.graphics.transformations import rotation_matrix
from graphics import MolecularGraphics

from interactionmanager import InteractionManager
from signalmanager import OneDSignalManager, TwoDSignalManager


class ChemLabMinimiser:
    def __init__(self):
        self.interaction_manager = get_interaction_manager(*get_twod_signal_manager().get_interaction_data())
        self.best_response_value = 1000000.0
        initial_coordinates = np.random.rand(self.interaction_manager.get_number_atoms(), 3)
        self.coordinate_manager = CoordinateManager(self, initial_coordinates, self.interaction_manager)
        self.iteration_number = 0
        self.fragments = self.generate_fragments_list()
        self.graphics = MolecularGraphics(self.coordinate_manager.get_coordinates(),
                                          self.interaction_manager.get_type_array(),
                                          np.array(self.interaction_manager.get_bonds())
                                          )
        self.hmbc_complete = False

    #TODO: Merge generate functions together
    def generate_fragments_list(self):
        """
        Creates a list of fragment objects
        """
        fragments = []
        unselected = [x for x in range(self.interaction_manager.get_number_atoms())]
        selected = []
        while len(unselected) > 0:
            new_atoms = set()
            new_atoms.add(unselected[0])
            global_indices = []
            global_bond_indices = []
            while len(new_atoms) > 0:
                a = list(new_atoms)[0]
                global_indices.append(a)
                new_atoms.remove(a)
                for bond_index, bond in enumerate(self.interaction_manager.get_bonds()):
                    if global_indices[-1] == bond[0] and bond[1] not in global_indices:
                        new_atoms.add(bond[1])
                    if global_indices[-1] == bond[1] and bond[0] not in global_indices:
                        new_atoms.add(bond[0])
                    if a in bond and list(bond) not in global_bond_indices:
                        global_bond_indices.append(list(bond))
            global_indices.sort()
            fragment = Fragment(global_indices, global_bond_indices, self.coordinate_manager, self.interaction_manager)
            fragments.append(fragment)
            selected += global_indices
            unselected = [x for x in range(self.interaction_manager.get_number_atoms()) if x not in selected]
        return fragments


    def main(self):
        mass_input = None
        while mass_input is None:
            try:
                mass_input = float(input("What is the Mr of this molecule: "))
            except ValueError:
                print("Mr Must be a floating point number")
                continue

        mass_hydrocarbon = self.interaction_manager.get_total_mass()
        print(self.coordinate_manager.calculate_response_value(True))
        self.interaction_manager.print_matrix()
        print("Total mass is %s, current mass is %s, delta mass is %s" % (mass_input, mass_hydrocarbon, mass_input-mass_hydrocarbon))
        print(self.coordinate_manager.calculate_response_value())
        self.graphics.run(self.iteration)

    def infer_hmbc(self):
        hmbc_interactions = self.interaction_manager.get_all_interactions(interaction_type=3)
        hmbc_interactions = [x for x in hmbc_interactions if self.bond_distance(x[0], x[1]) not in [2, 3]]

        if len(hmbc_interactions) == 0:
            print("###\n"*5 + "All HMBC interations satisfied\n"+"###\n"*5)
            self.hmbc_complete = True
            return
        else:
            print("Unsatisfied HMBCs: %s"% str(hmbc_interactions))

            """
            a = "~"
            while a != "":
                try:
                    return
                    a = input(">>>")
                    eval(a)
                except:
                    break
            """

        for hmbc in hmbc_interactions:
            for frag1 in self.fragments:
                for frag2 in self.fragments:
                    if frag1 != frag2 or len(self.fragments) == 1:
                        hw = min(hmbc)
                        cz = max(hmbc)
                        if not(hw in frag1.global_indices and cz in frag2.global_indices) and not(cz in frag1.global_indices and hw in frag2.global_indices):
                            continue
                        if self.interaction_manager.get_type(hw) != "H":
                            raise Exception("HMBC interaction invalid! Should he H-C interaction")

                        bonded_carbon_bonds = [x for x in self.interaction_manager.get_bonds() if hw in x]
                        if len(bonded_carbon_bonds) > 1:
                            raise Exception("Something Wrong Here")

                        if len(bonded_carbon_bonds) == 0:
                            print("%s has =/= 1 number of bonded carbons: %s" % (hw, str(bonded_carbon_bonds)))
                            continue

                        if len(bonded_carbon_bonds) == 1:
                            cx = max(bonded_carbon_bonds[0])

                        bond_formed = self.add_hmbc_bond(cx, cz, frag1, frag2, hmbc) or self.add_hmbc_bond(cz, cx, frag2, frag1, hmbc)
                        if bond_formed:
                            return

    def bond_distance(self, i, j):
        fragment = [x for x in self.fragments if i in x.global_indices][0]

        one_bond_atoms = self.get_adjacent(fragment, i)
        two_bond_atoms = uniq([self.get_adjacent(fragment, x, one_bond_atoms+[i]) for x in one_bond_atoms])
        three_bond_atoms = uniq([self.get_adjacent(fragment, x, one_bond_atoms+two_bond_atoms+[i]) for x in two_bond_atoms])

        if i == j:
            return 0
        if j in one_bond_atoms:
            return 1
        if j in two_bond_atoms:
            return 2
        if j in three_bond_atoms:
            return 3


    def get_adjacent(self, fragment, i, expanded=[]):
        adjacent_bonds = [x for x in self.interaction_manager.get_bonds() if i in x]
        adjacent_atoms = [x for x in uniq(adjacent_bonds) if x != i and x not in expanded]
        return adjacent_atoms

    def add_hmbc_bond(self, cx, cz, frag1, frag2, hmbc):
        free_valencies = self.interaction_manager.get_free_valencies()
        if free_valencies[cz] == 0:
            return False

        if free_valencies[cx] == 0:
            cx_adjacent_carbons = [list(x) for x in self.interaction_manager.get_bonds() if cx in x]
            cx_adjacent_carbons = [x for x in uniq(cx_adjacent_carbons) if x!=cx and self.interaction_manager.get_free_valencies()[x]!=0]
            print("Partially valent carbons adjacent to %s: %s" % (cx, str(cx_adjacent_carbons)))
            if len(cx_adjacent_carbons) == 1:
                cy = cx_adjacent_carbons[0]
                bond = [cy, cz]
                if cy == cz or bond in self.interaction_manager.get_bonds():
                    return False
                print("HMBC:", hmbc)
                print("Bond Chain: ", cx, cy, cz)
                self.add_bond(frag1, frag2, bond)
                return True
            else:
                return False

    def iteration(self):
        for fragment in self.fragments:
            fragment.project_bond_lengths()

        print("Starting Iteration %s with %s fragments"%(self.iteration_number, len(self.fragments)))
        if not self.hmbc_complete:
            self.infer_hmbc()

        self.iteration_number += 1

        for fragment in self.fragments:
            fragment.start_iter()
            fragment.project_bond_lengths()
            fragment.rotate_bonds(sequence=1)
            fragment.rotate_bonds(sequence=1)
            fragment.verify()
            fragment.start_iter()
            fragment.translate()
            fragment.verify()
            self.update_coordinates()

        self.graphics.update_coordinates(self.coordinate_manager.get_coordinates())
        self.interaction_manager.write_numpy_to_mol("tempfile.mol", self.coordinate_manager.get_coordinates())




    def add_bond(self, fragment1, fragment2, bond):
        self.interaction_manager.set_interaction(bond[0], bond[1], 4)
        print("Adding Bond %s"%str(bond))

        self.interaction_manager.add_bond(bond, 1)
        if len(self.fragments) > 1:
            self.merge_fragments(fragment1, fragment2, bond)
            self.fragments.remove(fragment1)
            if fragment1 != fragment2:
                self.fragments.remove(fragment2)

        self.coordinate_manager.reset()
        self.update_graphics()

    def update_graphics(self):
        self.graphics.update(self.coordinate_manager.get_coordinates(), self.interaction_manager.get_type_array(), self.interaction_manager.get_bonds())

    def update_coordinates(self):
        self.graphics.update_coordinates(self.coordinate_manager.get_coordinates())

    def merge_fragments(self, fragment1, fragment2, new_bond):
        new_fragment_indices = fragment1.global_indices + fragment2.global_indices
        new_fragment_bonds = fragment1.global_bond_indices + fragment2.global_bond_indices+[new_bond]
        merged_fragment = Fragment(new_fragment_indices, new_fragment_bonds, self.coordinate_manager, self.interaction_manager)
        self.fragments.append(merged_fragment)
        return merged_fragment


class Fragment:
    def __init__(self, global_indices, global_bond_indices, coordinate_manager, interaction_manager):
        self.global_indices = global_indices
        self.global_bond_indices = global_bond_indices
        self.coordinate_manager = coordinate_manager
        self.interaction_manager = interaction_manager
        self.project_bond_lengths()
        self.best_response_value = 1000
        self.verify_index = 0
        self.has_hydrogens = [False for x in range(len(self.coordinate_manager.atom_coordinates))]
        for x in self.global_bond_indices:
            if self.interaction_manager.get_type(x[0]) == "H":
                self.has_hydrogens[x[1]] = True
            if self.interaction_manager.get_type(x[1]) == "H":
                self.has_hydrogens[x[0]] = True

    def project_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.global_bond_indices:

            if self.interaction_manager.get_interaction(*bond) == 9 or bond[0] == bond[1]:
                continue
            frozen, unfrozen = self.bisect_on_bond(bond)
            bond_vector = atom_coordinates[bond[1]] - atom_coordinates[bond[0]]
            current_bond_length = np.linalg.norm(bond_vector)
            bond_length = self.interaction_manager.get_bond_length(*bond)
            if abs(current_bond_length-bond_length) > 0.01:
                corrected_bond_vector = bond_vector * (bond_length/current_bond_length)
                atom_coordinates[unfrozen] += (corrected_bond_vector-bond_vector)
        self.coordinate_manager.update(atom_coordinates)

    def reduce_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.global_bond_indices:
            if self.interaction_manager.get_interaction(*bond)==9:
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
        unfrozen_atoms = [x for x in uniq(bonds) if x not in frozen_atoms]#


        return list(frozen_atoms), unfrozen_atoms

    def rotate_bonds(self, bond_index=None, angle=None, sequence=1):
        #self.reduce_bond_lengths()
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

            if random.randrange(0, 3) == 0:
                random_vector = np.random.rand(3)
                rotation_axis = np.cross(point_b-point_a, random_vector)
            elif random.randrange(0, 3) == 1:
                rotation_axis = np.random.rand(3)
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
        self.temperature = 10

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

    carbons = len([signal for signal in oned_signal_manager.signals if signal.signal_type == "C"])
    hydrogens = len([signal for signal in oned_signal_manager.signals if signal.signal_type == "H"])
    print("%s Carbons and %s Hydrogens"%(carbons, hydrogens))

    input("Press Enter to Continue...")


    twod_signal_manager = TwoDSignalManager(oned_signal_manager)
    twod_signal_manager.add_nmr_signals("COSY", 'resources/nmr/twod/cosy/cosy_peak_data.txt')
    twod_signal_manager.add_nmr_signals("HSQC", 'resources/nmr/twod/hsqc/hsqc_peak_data.txt')
    twod_signal_manager.add_nmr_signals("HMBC", 'resources/nmr/twod/hmbc/hmbc_peak_data.txt')
    return twod_signal_manager

def get_interaction_manager(interaction_matrix, type_array, shift_data):
    interaction_manager = InteractionManager(1001, interaction_matrix, type_array, shift_data)
    interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=0.2, depth=10)
    interaction_manager.add_new_interaction(index=2, interaction_name="HSQC 1-Bond H-C          ", repulsive_amplitude=0.8, repulsive_time_constant=2.26, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_hmbc_interaction(index=3, interaction_name="HMBC 2/3 Bond H-C        ", repulsive_amplitude=0.8, repulsive_time_constant=6.50, depth=3.8, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=4, interaction_name="INAD Single 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.49, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=5, interaction_name="INAD Double 1-Bond C=C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.15, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=6, interaction_name="NOESY Normalised Function", repulsive_amplitude=0.8, repulsive_time_constant=2.03, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    return interaction_manager


def uniq(x):
    x = [list(y) for y in x]
    return list(set(sum(x, [])))

def end():
    input("Press Any Key to Quit")
    raise SystemExit

def main():
    clm = ChemLabMinimiser()
    clm.main()
    end()

if __name__ == "__main__":
    #import cProfile
    #cProfile.run('main()')
    main()
