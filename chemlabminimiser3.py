import numpy as np
import random
from chemlab.graphics.transformations import rotation_matrix
import os

from graphics import MolecularGraphics
from interactionmanager import InteractionManager
from signalmanager import get_twod_signal_manager, InteractionValues



def clear():
    os.system('cls' if os.name == 'nt' else 'clear')

class ChemLabMinimiser:
    def __init__(self, path):
        self.interaction_manager = get_interaction_manager(path)
        #self.interaction_manager.plot_all_interactions()
        self.best_response_value = 1000000.0
        initial_coordinates = np.random.rand(self.interaction_manager.get_number_atoms(), 3)
        self.coordinate_manager = CoordinateManager(initial_coordinates)
        self.iteration_number = 0
        self.fragments = self.generate_fragments_list()
        self.graphics = MolecularGraphics(self.coordinate_manager.get_coordinates(),
                                          self.interaction_manager
                                          )
        self.hmbc_complete = False

        def viewer_addbond(atom1=None, atom2=None, bond_order=1, interaction_type=0):
            if not atom1 or not atom2:
                selected_atoms = self.get_selected_atoms()
                if len(selected_atoms) != 2:
                    print("Select Exactly Two Atoms to Bond")
                    return
                atom1, atom2 = selected_atoms[0], selected_atoms[1]
            self.add_bond(atom1, atom2)
            self.fragments = self.generate_fragments_list()
            self.update_graphics()
        self.graphics.add_bondsetter(viewer_addbond)

        self.graphics.add_atomsetter(self.add_atom)





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

    def check_unsatisfied_hmbc(self):
        hmbc_interactions = self.interaction_manager.get_all_interaction_atoms(InteractionValues.HMBC)

        for hmbc in hmbc_interactions:
            if self.bond_distance(hmbc[0], hmbc[1]) in [2, 3]:
                self.interaction_manager.set_interaction(hmbc[0].index_value, hmbc[1].index_value, InteractionValues.DEFAULT)
                hmbc_interactions.remove(hmbc)
                self.update_graphics()
            elif self.bond_distance(hmbc[0], hmbc[1]) == 4:
                raise Exception("4-Bond HMBC between signals at %s %s" % (hmbc[0].shift_value, hmbc[1].shift_value))
        return hmbc_interactions

    def add_atom(self, atom_type):

        self.interaction_manager.add_atom(0, atom_type)
        self.coordinate_manager.add_atom()
        self.fragments = self.generate_fragments_list()
        self.update_graphics()

    def infer_hmbc(self):

        hmbc_interactions = self.check_unsatisfied_hmbc()

        if len(hmbc_interactions) == 0:
            print("All HMBC interations satisfied")
            return
        else:
            print("Unsatisfied HMBCs: %s" % str(hmbc_interactions))
            if len(hmbc_interactions) == 2:
                return


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
        four_bond_atoms = uniq([self.get_adjacent(x, one_bond_atoms+two_bond_atoms+three_bond_atoms+[atom1]) for x in three_bond_atoms])

        if atom1 == atom2:
            return 0
        if atom2 in one_bond_atoms:
            return 1
        if atom2 in two_bond_atoms:
            return 2
        if atom2 in three_bond_atoms:
            return 3
        if atom2 in four_bond_atoms:
            return 4

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
                    self.add_bond(cy, cz, 1, InteractionValues.HMBC)

                    return True
            else:
                return False

    def get_selected_atoms(self):
        return self.graphics.get_selected_atoms()

    def iteration(self):
        clear()
        print("Hydrogens",self.interaction_manager.get_type_array().count("H"))
        print("Carbons",self.interaction_manager.get_type_array().count("C"))
        self.iteration_number += 1
        print("Starting Iteration %s with %s fragments" % (self.iteration_number, len(self.fragments)))
        selected_atoms = self.get_selected_atoms()

        for atom in selected_atoms:
            atom.print_data()

        self.infer_hmbc()
        self.interaction_manager.update_all_bonds()

        for fragment in self.fragments:
            fragment.project_bond_lengths()

        for fragment in self.fragments:
            self.start_iter()
            fragment.rotate_bonds()
            self.verify()

            self.start_iter()
            fragment.translate()
            self.verify()

            self.update_coordinates()

        self.interaction_manager.write_numpy_to_mol("tempfile.mol", self.coordinate_manager.get_coordinates())




    def add_bond(self, atom1, atom2, bond_order=1, interaction_type=0):
        self.best_response_value = 1000000
        self.interaction_manager.add_bond(atom1, atom2, bond_order, interaction_type)
        self.fragments = self.generate_fragments_list()
        print("Adding Bond %s" % str([atom1, atom2]))


    def update_graphics(self):
        self.graphics.update(self.coordinate_manager.get_coordinates())

    def update_coordinates(self):
        self.graphics.update_coordinates(self.coordinate_manager.get_coordinates())

    def main(self):
        mass_input = 0
        """
        mass_input = None
        while mass_input is None:
            try:
                mass_input = float(input("What is the Mr of this molecule: "))
            except ValueError:
                print("Mr Must be a floating point number")
        """

        current_molecular_mass = self.interaction_manager.get_total_mass()
        delta_mass = mass_input - current_molecular_mass

        print("Total mass is %s, current mass is %s, delta mass is %s" %
              (mass_input, current_molecular_mass, delta_mass))

        print(self.calculate_response_value())
        self.interaction_manager.print_matrix()
        self.graphics.run(self.iteration)

    def verify(self):

        response = self.calculate_response_value()
        if response < self.best_response_value:
            self.best_response_value = response
            print("New Best Response: %s "%(response))
            return True
        else:
            self.coordinate_manager.revert()
            return False

    def calculate_response_value(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        response = 0
        for i, v1 in enumerate(atom_coordinates):
            for j, v2 in enumerate(atom_coordinates):
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, v1, v2)
        return response

    def start_iter(self):
        self.coordinate_manager.start_iter()



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
                corrected_bond_vector = bond_vector * (bond_length/current_bond_length)
                atom_coordinates[frozen] += (corrected_bond_vector-bond_vector)

        self.coordinate_manager.update(atom_coordinates)

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



class CoordinateManager:
    def __init__(self, atom_coordinates):
        self.atom_coordinates = atom_coordinates
        self.old_atom_coordinates = np.copy(atom_coordinates)

    def start_iter(self):
        self.old_atom_coordinates = np.copy(self.atom_coordinates)

    def update(self, atom_coordinates):
        atom_coordinates -= atom_coordinates[-1]
        self.atom_coordinates = atom_coordinates

    def revert(self):
        self.atom_coordinates = self.old_atom_coordinates

    def get_coordinates(self):
        return np.copy(self.atom_coordinates)

    def add_atom(self):
        self.atom_coordinates = np.insert(self.atom_coordinates, len(self.atom_coordinates), values=0.1, axis=0)
        self.old_atom_coordinates = np.copy(self.atom_coordinates)


def get_interaction_manager(path):
    twod_signal_manager = get_twod_signal_manager(path)
    interaction_matrix, type_array, shift_data = twod_signal_manager.get_interaction_data()
    interaction_manager = InteractionManager(1001, interaction_matrix, type_array, shift_data)
    interaction_manager.add_default_interaction(InteractionValues.DEFAULT, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=0.2, depth=100)
    interaction_manager.add_hmbc_interaction(index=InteractionValues.HMBC,    interaction_name="HMBC 2/3 Bond H-C        ", repulsive_amplitude=0.8, repulsive_time_constant=6.50, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_spatial_interaction(index=InteractionValues.NOESY, interaction_name="NOESY Normalised Function", repulsive_amplitude=0.8, repulsive_time_constant=2.03, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    return interaction_manager

def uniq(x):
    x = [list(y) for y in x]
    return list(set(sum(x, [])))

def end():
    raise SystemExit

def main():
    path = "C:/users/martin/desktop/nmr samples/edmpc_new_format.zip"
    clm = ChemLabMinimiser(path)
    clm.main()
    end()

if __name__ == "__main__":
    main()
