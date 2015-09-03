import numpy as np
import random
from chemlab.graphics.transformations import rotation_matrix
import os
from graphics import MolecularGraphics
from interactionmanager import InteractionManager
from signalmanager import get_twod_signal_manager, InteractionValues
import time
from tkinter import filedialog


def clear():
    os.system('cls' if os.name == 'nt' else 'clear')


class EnergyMinimiser:
    def __init__(self, path=None, mass=None):
        """
        Constructor for the Minimiser Class
        Args:
            path (str) : Path to zip folder containing NMR files [Optional]
            mass (float) : Floating Point Value of Molecular Mass [Optional]
        Returns:
            None
        """
        if path is None:
            path = filedialog.askopenfilename(title="Open Zipped NMR File")
        self.interaction_manager = get_interaction_manager(path)

        self.mass_input = mass
        while self.mass_input is None:
            try:
                self.mass_input = float(input("What is the Mr of this molecule: "))
            except ValueError:
                print("Mr Must be a floating point number")

        self.best_response_value = 1000000.0
        self.iteration_number = 0
        self.frame_number = 0
        self.fragments = []
        self.fragments = self.generate_fragments_list()
        self.hmbc_complete = False
        self.paused = False
        self.update = False
        self.start_time = None
        self.graphics = MolecularGraphics(self.interaction_manager)
        clear()

        def viewer_add_bond(atom1=None, atom2=None, bond_order=1):
            """
            Function for the iPython terminal to manually add new bonds between atoms
            Can accept manually defined atoms, but usually used by selecting two atoms
            using the graphics interface and using this bond function interactively.

            Args:
                atom1 (Atom) : First Atom in the new bond
                atom2 (Atom) : Second Atom in the new bond
                bond_order (float): Order of the new bond
            Returns:
                None
            """
            if not atom1 or not atom2:
                selected_atoms = self.get_selected_atoms()
                if len(selected_atoms) != 2:
                    print("Select Exactly Two Atoms to Bond")
                    return
                atom1, atom2 = selected_atoms[0], selected_atoms[1]

            message = "Added Bond Between atoms %s and %s with bond order %s" % \
                      (atom1.index_value, atom2.index_value, bond_order)
            self.add_bond(atom1, atom2)
            self.fragments = self.generate_fragments_list()
            self.update = True
            unpause()
            save_state(message)
            for x in range(len(self.fragments)+1):
                self.iteration()
            pause()

        def viewer_add_atom(atom_type):
            """
            Function for the iPython terminal to manually add new atoms.
            Requires the new atom type to be defined. The Atom is added to the origin,
            and if another atom is selected, then a new bond is created between this
            atoms and the new atom

            Args:
                atom_type (Str): Type of Atom "C", "H", "N" etc
            Returns:
                Nones
            """
            selected_atoms = self.get_selected_atoms()
            message = "Added Atom of Type %s" % atom_type
            self.add_atom(atom_type)
            if len(selected_atoms) == 1:
                viewer_add_bond(selected_atoms[0], self.interaction_manager.atoms[-1])
            save_state(message)
            for x in range(len(self.fragments)+1):
                self.iteration()
            pause()

        def save_state(message="Manual State Save"):
            self.interaction_manager.savestate(message)

        def load_stage(state_number=None, revert=False):
            if state_number:
                state_number -= 1
            pause()
            self.best_response_value = 1000
            self.interaction_manager.loadstate(state_number, revert)
            self.update = True
            self.fragments = self.generate_fragments_list()
            unpause()

        def revert():
            load_stage(0, True)

        def pause():
            self.paused = True

        def unpause():
            self.frame_number = 0
            self.start_time = time.time()
            self.paused = False

        """Injecting Interactive Functions into IPython Terminal"""
        namespace = self.graphics.get_namespace()
        namespace.add_bond = viewer_add_bond
        namespace.add_atom = viewer_add_atom
        namespace.save_state = save_state
        namespace.revert = revert
        namespace.load_state = load_stage
        namespace.pause = pause
        namespace.unpause = unpause
        namespace.resume = unpause

    def generate_fragments_list(self):
        """
        Traverses the bonds between atoms in the system and
        generates fragment objects of contiguous chains of atoms

        Args:
            None
        Returns:
            fragments (List[Fragment]) : List of Fragment objects

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
            new_fragment = Fragment(new_fragment_atoms, self.interaction_manager)
            fragments.append(new_fragment)
            selected += new_fragment_atoms
            unselected = [x for x in self.interaction_manager.atoms if x not in selected]
        return fragments

    def check_unsatisfied_hmbc(self):
        """
        Gets all HMBC interactions from the interaction manager and checks which are satisfied

        Args:
            None
        Returns:
            hmbc_interactions ( List[[x,y]] ) : List of unsatisfied HMBC interactions defined by atom indices
        """
        hmbc_interactions = self.interaction_manager.get_all_interaction_atoms(InteractionValues.HMBC)
        for hmbc in hmbc_interactions:
            if self.bond_distance(hmbc[0], hmbc[1]) in [2, 3]:
                self.interaction_manager.set_interaction(hmbc[0].index_value, hmbc[1].index_value, InteractionValues.DEFAULT)
                hmbc_interactions.remove(hmbc)
                self.update = True
            elif self.bond_distance(hmbc[0], hmbc[1]) == 4:
                raise Exception("4-Bond HMBC between signals at %s %s" % (hmbc[0].shift_value, hmbc[1].shift_value))
        return hmbc_interactions

    def add_atom(self, atom_type):
        """
        Manually adds a new atom to the system, the regenerates the fragments
        Fragment list is regenerated to incorporate the new atom

        Args:
            atom_type (str) : Type of Atom (e.g. "H", "C", "N", etc

        """
        self.interaction_manager.add_atom(0, atom_type)
        self.fragments = self.generate_fragments_list()
        self.update = True

    def infer_hmbc(self):
        """
        Logically determines which bonds can be inferred with certainty
        according to the 2/3 Bond HMBC interactions, by considering those
        which cannot be 2-Bond interactions and finding the only
        configuration that would satisfy a 3-Bond correlation.

        Args:
            None
        Returns:
            None
        """

        self.hmbc_complete = False
        hmbc_interactions = self.check_unsatisfied_hmbc()
        if len(hmbc_interactions) == 0:
            self.hmbc_complete = True

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
                            raise Exception("Error: Hydrogen bonded to multiple carbons")
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
        """
        Finds the bond distance between adjacent atoms in the system.
        Whilst this can be generalised to n-bonds, we only care about bonds in range 1-4
         => Stops counting bonds beyond this due to unnecessary computation

        Args:
            atom1 (Atom) : First Atom
            atom2 (Atom) : Second Atom
        Returns:
            distance (int) : Distance between atoms in bonds

        """
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

    def printout(self):
        clear()
        end_section = "|____________________|\n\n"

        print("Attempting to solve structure" + "." * (self.iteration_number % 4))
        print("\n")


        print("Starting Iteration %s with %s fragments" % (self.iteration_number, len(self.fragments)))
        print("\n")

        if self.mass_input > 0:
            current_molecular_mass = self.interaction_manager.get_total_mass()
            delta_mass = self.mass_input - current_molecular_mass
            print("___Mass Information___")
            print("Total mass: %s" % self.mass_input)
            print("Current mass: %s" % current_molecular_mass)
            print("Delta mass is: %s" % delta_mass)
            print(end_section)

        selected_atoms = self.get_selected_atoms()
        print("____Selected Atoms____")
        for atom in selected_atoms:
            atom.print_data()
        print(end_section)

        print("______Atom Types______")
        atom_types = self.interaction_manager.get_type_array()
        for atom_type in set(atom_types):
            print("Number of %s: " % atom_type, atom_types.count(atom_type))
        print(end_section)

        print("____Response Value____")
        print("Current best Response: %s " % self.best_response_value)
        print(end_section)

        if not self.paused:
            print("___Iteration Speed____")
            running_time = time.time() - self.start_time
            iteration_rate = round(self.frame_number/running_time,2)
            print("Current Iteration Rate: %s " % iteration_rate)
            print(end_section)

        print("__HMBC Interactions__")
        if self.hmbc_complete:
            print("All HMBC Interactions Satisfied")
        else:
            hmbcs = self.interaction_manager.get_all_interaction_atoms(InteractionValues.HMBC)
            print("Unsatisfied HMBCS:")
            print(hmbcs)
        print(end_section)

    def iteration(self):
        try:
            if self.paused:
                self.printout()
                return
            self.printout()
            self.iteration_number += 1
            self.frame_number += 1

            for x in range(0, 5):
                self.infer_hmbc()
                self.interaction_manager.update_all_bonds()

                fragment = self.fragments[0]
                if True:
                    fragment.project_bond_lengths()
                    self.start_iter()
                    fragment.rotate_bonds()
                    self.verify()

                    self.start_iter()
                    fragment.translate()
                    self.verify()
                if self.update and self.iteration_number > 25:
                    self.update_graphics()
                    self.update = False
                else:
                    self.update_coordinates()
                self.fragments = self.fragments[1:] + [fragment]
                self.interaction_manager.write_numpy_to_mol("resources/tempfile.mol")
        except Exception as e:
            print(str(e))
            self.update_graphics()

    def start_iter(self):
        self.interaction_manager.start_iter()

    def verify(self):
        response = self.interaction_manager.calculate_response_value()
        if response < self.best_response_value:
            self.best_response_value = response
            return True
        else:
            self.interaction_manager.revert()
            return False

    def add_bond(self, atom1, atom2, bond_order=1, interaction_type=0):
        self.best_response_value = 1000000
        self.interaction_manager.add_bond(atom1, atom2, bond_order, interaction_type)
        self.fragments = self.generate_fragments_list()
        print("Adding Bond %s" % str([atom1, atom2]))
        self.update = True

    def update_graphics(self):
        print("Updating Graphics")
        self.frame_number = 1
        self.start_time = time.time()
        self.graphics.update()

    def update_coordinates(self):
        self.graphics.update_coordinates()

    def main(self):
        self.start_time = time.time()
        self.graphics.run(self.iteration)

    def calculate_response_value(self):
        return self.interaction_manager.calculate_response_value()


class Fragment:
    def __init__(self, atoms, interaction_manager):
        self.atoms = list(atoms)
        self.bonds = uniq([atom.bonds for atom in self.atoms])
        self.coordinate_manager = interaction_manager.coordinate_manager
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
            corrected_bond_vector = bond_vector * (bond_length/current_bond_length)
            atom_coordinates[frozen] += (corrected_bond_vector-bond_vector)
        self.coordinate_manager.update(atom_coordinates)

    def translate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        first_atom = new_atom_coordinates[self.atoms[0].index_value]
        for x in [0, 1, 2]:
            if first_atom[x] > 1:
                delta = first_atom[x] - 1
                vec = [0, 0, 0]
                vec[x] = delta
                new_atom_coordinates[[y.index_value for y in self.atoms]] -= vec



        xvec = np.random.normal(loc=1.0, scale=3.0)*(np.random.rand(3) - 0.5)
        xvec /= np.linalg.norm(xvec)
        vec = random.uniform(0, 0.6)*xvec
        new_atom_coordinates[[y.index_value for y in self.atoms]] += vec
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
            for bond in self.bonds:

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


def get_interaction_manager(path):
    twod_signal_manager = get_twod_signal_manager(path)
    interaction_matrix, type_array, shift_data, distance_matrix = twod_signal_manager.get_interaction_data()
    interaction_manager = InteractionManager(1001, interaction_matrix, type_array, shift_data, distance_matrix)

    interaction_manager.add_default_interaction(index=InteractionValues.DEFAULT,
                                                interaction_name="Default Repulsive",
                                                repulsive_amplitude=0.8,
                                                repulsive_time_constant=0.2,
                                                depth=1000)

    interaction_manager.add_hmbc_interaction(index=InteractionValues.HMBC,
                                             interaction_name="HMBC 2/3 Bond H-C",
                                             repulsive_amplitude=0.8,
                                             repulsive_time_constant=6.50,
                                             depth=3,
                                             attractive_amplitude=0.6,
                                             attractive_time_constant=200,
                                             power=3)

    interaction_manager.add_spatial_interaction(index=InteractionValues.NOESY,
                                                interaction_name="NOESY Normalised Function",
                                                repulsive_amplitude=0.8,
                                                repulsive_time_constant=2.03,
                                                depth=3,
                                                attractive_amplitude=0.6,
                                                attractive_time_constant=200,
                                                power=3)
    return interaction_manager


def uniq(x):
    x = [list(y) for y in x]
    return list(set(sum(x, [])))


def main():
    import sys
    mass = path = None
    if len(sys.argv) > 1:
        path = sys.argv[1]
    if len(sys.argv) > 2:
        mass = int(sys.argv[2])
    clm = EnergyMinimiser(path=path, mass=mass)
    clm.main()


if __name__ == "__main__":
    main()


