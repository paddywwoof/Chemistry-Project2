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
                                          self.interaction_manager
                                          )
        self.hmbc_complete = False

    #TODO: Merge generate functions together
    def generate_fragments_list(self):
        """
        Creates a list of fragment objects
        """
        fragments = []
        unselected = self.interaction_manager.atoms
        selected = []

        while len(unselected) > 0:

            fringe_atoms = set()
            fringe_atoms.add(unselected[0])
            new_fragment_atoms = set()
            while len(fringe_atoms) > 0:
                new_atom = list(fringe_atoms)[0]
                new_fragment_atoms.add(new_atom)
                fringe_atoms.remove(new_atom)
                for atom in new_atom.get_adjacent():
                    if atom not in new_fragment_atoms:
                        fringe_atoms.add(atom)
            new_fragment = Fragment(new_fragment_atoms, self.coordinate_manager, self.interaction_manager)
            fragments.append(new_fragment)
            selected += new_fragment_atoms
            unselected = [x for x in self.interaction_manager.atoms if x not in selected]
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
        hmbc_interactions = self.interaction_manager.get_all_interaction_atoms(interaction_type=3)

        graphics_update = False
        for hmbc in hmbc_interactions:
            if self.bond_distance(hmbc[0], hmbc[1]) in [2, 3]:
                graphics_update = True
                self.interaction_manager.set_interaction(hmbc[0].index_value, hmbc[1].index_value, 0)
                hmbc_interactions.remove(hmbc)
        if graphics_update:
            self.update_graphics()

        #hmbc_interactions = [x for x in hmbc_interactions if self.bond_distance(x[0], x[1]) not in [2, 3]]

        if len(hmbc_interactions) == 0:
            print("###\n"*5 + "All HMBC interations satisfied\n"+"###\n"*5)
            self.hmbc_complete = True
            return
        else:
            print("Unsatisfied HMBCs: %s" % str(hmbc_interactions))
            print("Free Valencies", self.interaction_manager.get_free_valencies())

        for hmbc in hmbc_interactions:
            for frag1 in self.fragments:
                for frag2 in self.fragments:
                    if frag1 != frag2 or len(self.fragments) == 1:
                        if hmbc[0].atom_type == "H" and hmbc[1].atom_type == "C":
                            hw, cz = hmbc
                        elif hmbc[1].atom_type == "H" and hmbc[0].atom_type == "C":
                            cz, hw = hmbc
                        else:
                            raise Exception("HMBC should be H-C interaction")
                        if not(hw in frag1.atoms and cz in frag2.atoms) and not(hw in frag2.atoms and cz in frag1.atoms):
                            continue
                        bonded_carbons = hw.get_adjacent()

                        if len(bonded_carbons) > 1:
                            raise Exception("Something Wrong Here")
                        elif len(bonded_carbons) == 0:
                            continue
                        else:
                            cx = bonded_carbons[0]

                        hmbc_type1 = self.add_hmbc_bond(cx, cz, frag1, frag2, hmbc)
                        if hmbc_type1:
                            return True

                        hmbc_type2 = self.add_hmbc_bond(cz, cx, frag2, frag1, hmbc)
                        if hmbc_type2:
                            return True

    def bond_distance(self, atom1, atom2):
        one_bond_atoms = self.get_adjacent(atom1)
        two_bond_atoms = uniq([self.get_adjacent(x, one_bond_atoms+[atom1]) for x in one_bond_atoms])
        three_bond_atoms = uniq([self.get_adjacent(x, one_bond_atoms+two_bond_atoms+[atom1]) for x in two_bond_atoms])

        if atom1 == atom2:
            return 0
        if atom2 in one_bond_atoms:
            return 1
        if atom2 in two_bond_atoms:
            return 2
        if atom2 in three_bond_atoms:
            return 3


    def get_adjacent(self, atom, expanded=[]):
        return [adjacent_atom for adjacent_atom in atom.get_adjacent() if adjacent_atom not in expanded]

    def add_hmbc_bond(self, cx, cz, frag1, frag2, hmbc):
        if cz.get_free_valency() == 0:
            return False
        if cx.get_free_valency() == 0:
            cx_adjacent_carbons = [x for x in cx.get_adjacent() if x.get_free_valency() != 0]
            if len(cx_adjacent_carbons) == 1:
                cy = cx_adjacent_carbons[0]
                if cy == cz or cy in cz.get_adjacent() or cz in cy.get_adjacent():
                    return False
                else:
                    print("HMBC:", hmbc)
                    print("Bond Chain: ", cx, cy, cz)
                    self.add_bond(cy, cz, 1, 3)
                    return True
            else:
                return False

    def iteration(self):
        print("Starting Iteration %s with %s fragments"%(self.iteration_number, len(self.fragments)))
        if not self.hmbc_complete:
            if self.infer_hmbc():
                self.update_graphics()
        self.interaction_manager.update_all_bonds()

        self.iteration_number += 1

        for fragment in self.fragments:
            fragment.project_bond_lengths()


        for fragment in self.fragments:
            if self.iteration_number < 150:
                fragment.start_iter()
                fragment.rotate_bonds()
                fragment.rotate_bonds()
                fragment.verify()

                fragment.start_iter()
                fragment.translate()
                fragment.verify()

                fragment.start_iter()
                fragment.rotate()
                fragment.verify()


            self.update_coordinates()
        self.interaction_manager.write_numpy_to_mol("tempfile.mol", self.coordinate_manager.get_coordinates())

    def add_bond(self, atom1, atom2, bond_order, interaction_type):
        self.interaction_manager.add_bond(atom1, atom2, bond_order, interaction_type)
        self.fragments = self.generate_fragments_list()
        print("Adding Bond %s" % str([atom1, atom2]))
        self.coordinate_manager.reset()
        self.update_graphics()

    def update_graphics(self):
        self.graphics.update(self.coordinate_manager.get_coordinates(), self.interaction_manager)

    def update_coordinates(self):
        self.graphics.update_coordinates(self.coordinate_manager.get_coordinates())


class Fragment:
    def __init__(self, atoms, coordinate_manager, interaction_manager):
        self.atoms = list(atoms)
        self.bonds = uniq([atom.bonds for atom in self.atoms])
        self.coordinate_manager = coordinate_manager
        self.interaction_manager = interaction_manager
        self.project_bond_lengths()
        self.best_response_value = 1000
        self.verify_index = 0

    def project_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.bonds:
            frozen, unfrozen = self.bisect_on_bond(bond)
            bond_vector = atom_coordinates[bond.atom1.index_value] - atom_coordinates[bond.atom2.index_value]
            current_bond_length = np.linalg.norm(bond_vector)
            bond_length = bond.bond_length

            if abs(current_bond_length-bond_length) > 0.1:

                print(bond_length/current_bond_length)

                corrected_bond_vector = bond_vector * (bond_length/current_bond_length)
                atom_coordinates[frozen] += (corrected_bond_vector-bond_vector)

        self.coordinate_manager.update(atom_coordinates)

    def bisect_on_bond(self, bond, freeze_left=True):
        frozen_atoms = set()
        fringe_atoms = set()

        if freeze_left:
            fringe_atoms.add(bond.atom1)
            block_atom = bond.atom2
        else:
            fringe_atoms.add(bond.atom2)
            block_atom = bond.atom1

        while len(fringe_atoms) > 0:
            new_atom = list(fringe_atoms)[0]
            for bond in new_atom.bonds:
                if block_atom not in bond and new_atom not in frozen_atoms:
                    fringe_atoms.add(bond.other(new_atom))
            frozen_atoms.add(new_atom)
            fringe_atoms.remove(new_atom)

        unfrozen_indices = [x.index_value for x in self.atoms if x not in frozen_atoms]
        frozen_indices = [x.index_value for x in list(frozen_atoms)]
        return frozen_indices, unfrozen_indices

    def rotate_bonds(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        if len(self.bonds) > 0:
            bond = random.choice(self.bonds)

            index1 = bond.atom1.index_value
            index2 = bond.atom2.index_value

            if random.randrange(0, 2) == 0:
                index1, index2 = index2, index1

            point_a = np.copy(new_atom_coordinates[index1])
            point_b = np.copy(new_atom_coordinates[index2])

            uf = self.bisect_on_bond(bond, True)[1]
            angle = (np.pi/100) * random.uniform(0, 25)
            if random.randrange(0, 2) == 0:
                random_vector = np.random.rand(3)
                rotation_axis = np.cross(point_b-point_a, random_vector)
            else:
                rotation_axis = point_b - point_a

            rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
            try:
                rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
            except Exception as e:
                input(str(e))

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
        new_atom_coordinates[[x.index_value for x in self.atoms]] += vec

        self.coordinate_manager.update(new_atom_coordinates)

    def rotate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        atom = random.choice(new_atom_coordinates)
        new_atom_coordinates -= atom
        angle = np.random.normal(np.pi*0.5, np.pi*0.2)
        rotation_axis = np.random.rand(3)
        rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
        indices = [x.index_value for x in self.atoms]
        new_atom_coordinates[indices] = np.dot(new_atom_coordinates[indices], rotation.T)
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

        atom_coordinates -= atom_coordinates[0]

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
    interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=0.2, depth=100)
    interaction_manager.add_spatial_interaction(index=1, interaction_name="COSY 3-Bond H-H          ", repulsive_amplitude=0.8, repulsive_time_constant=6.28, depth=0, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_bond_interaction(index=2,    interaction_name="HSQC 1-Bond H-C          ", bond_length=0.109)
    interaction_manager.add_hmbc_interaction(index=3,    interaction_name="HMBC 2/3 Bond H-C        ", repulsive_amplitude=0.8, repulsive_time_constant=6.50, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_bond_interaction(index=4,    interaction_name="INAD Single 1-Bond C-C   ", bond_length=0.154)
    interaction_manager.add_spatial_interaction(index=6, interaction_name="NOESY Normalised Function", repulsive_amplitude=0.8, repulsive_time_constant=2.03, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_bond_interaction(index=7,    interaction_name="Inferred C=O Bond Type   ", bond_length=0.142)
    return interaction_manager


def uniq(x):
    x = [list(y) for y in x]
    return list(set(sum(x, [])))

def end():
    raise SystemExit

def main():
    clm = ChemLabMinimiser()
    clm.main()
    end()

if __name__ == "__main__":
    #import cProfile
    #cProfile.run('main()')
    main()
