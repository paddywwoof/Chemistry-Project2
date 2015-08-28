import numpy as np
from chemlab.graphics.renderers import BallAndStickRenderer
from chemlab.graphics import QtViewer
from chemlab.graphics.transformations import rotation_matrix
from chemlab.io.handlers import MolIO
import random

from filemanager import FileManager

global_frag_distances = {(12, 6): 4.63,
                         (10, 6): 4.63,
                         (7, 8): 1.78,
                         (10, 8): 2.84,
                         (10, 7): 2.99,
                         (10, 11): 2.45,
                         (11, 9): 2.43,
                         (9, 12): 2.46,
                         (10, 9): 4.25,
                         (11, 12): 4.25
                  }

class ChemLabMinimiser:
    def __init__(self):
        self.file_manager = FileManager()
        system, self.atom_coordinates = self.read_coordinates()
        #  self.interaction_manager = get_interaction_manager(get_twod_signal_manager().get_interaction_matrix(), system.type_array.tolist(), np.zeros(shape=system.type_array.shape))
        self.interaction_manager = get_interaction_manager(*get_twod_signal_manager().get_interaction_data())
        self.viewer = QtViewer()
        self.best_response_value = 1000000.0
        self.best_atom_coordinates = np.copy(self.atom_coordinates)
        self.fragments = self.generate_fragments(system)
        self.atom_coordinates = self.fix_hydrogen_bond_lengths(system, self.atom_coordinates)




    def generate_fragments(self, system):
        fragments = []
        groups = self.split_groups(system)
        for group in groups:
            new_bonds_reduced = np.array([[group[0].index(x[0]), group[0].index(x[1])] for x in group[1]])
            fragment = Fragment(group[0], system.type_array, np.array(group[1]),new_bonds_reduced, group[2], self.viewer, self.atom_coordinates, self.best_response_value, self.interaction_manager)
            fragments.append(fragment)
        return fragments

    def read_coordinates(self):
        infile = open('../outfile.mol', 'rb')
        mol_io = MolIO(infile)
        system = mol_io.read('molecule')
        atom_coordinates = system.r_array * 1.818

        return system, atom_coordinates

    def fix_hydrogen_bond_lengths(self, system, atom_coordinates):
        for bond in system.bonds:
            bond_types = [system.type_array[bond[0]], system.type_array[bond[1]]]
            if bond_types == ["C", "C"]:
                print("CC Bond: ", np.linalg.norm(atom_coordinates[bond[0]] - atom_coordinates[bond[1]]))
            if "H" in bond_types and "C" in bond_types:
                if bond_types[0] == "H":
                    h, c = bond
                else:
                    c, h = bond
                point_a = np.copy(atom_coordinates[c])
                point_b = np.copy(atom_coordinates[h])

                vec = (point_b-point_a)
                vec /= np.linalg.norm(vec)
                vec *= 0.109
                atom_coordinates[h][:] = point_a + vec
        return atom_coordinates

    def update(self):
        for x in [0, 1, 2]:
            self.atom_coordinates[:, x] -= min(self.atom_coordinates[:, x])
        for fragment in self.fragments:
            self.atom_coordinates, self.best_response_value = fragment.update(self.atom_coordinates, self.best_response_value)

    def main(self):
        self.interaction_manager.print_matrix()
        #self.interaction_manager.plot_all_interactions()
        for fragment in self.fragments:
            self.atom_coordinates = fragment.project_bond_lengths(self.atom_coordinates)
        """
        for fragment in self.fragments:
            self.atom_coordinates = fragment.reduce(self.atom_coordinates)
        """
        self.viewer.schedule(self.update, 1)
        self.viewer.run()
        raise SystemExit

    def split_groups(self, system):
        groups = []
        unselected = [x for x in range(system.n_atoms)]
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

                for bond_index, bond in enumerate(system.bonds):
                    if atoms[-1] == bond[0] and bond[1] not in atoms:
                        new_atoms.add(bond[1])
                    if atoms[-1] == bond[1] and bond[0] not in atoms:
                        new_atoms.add(bond[0])
                    if a in bond and bond.tolist() not in new_bonds:
                        new_bonds.append(bond.tolist())
                        new_bond_orders.append(system.bond_orders[bond_index])
            atoms.sort()
            groups.append([atoms, new_bonds, new_bond_orders])
            unselected = [x for x in range(system.n_atoms) if x not in sum([y[0] for y in groups],[])]
        return groups


class Fragment:
    def __init__(self, global_indices, type_array, bonds, bonds_reduced, bond_orders, v, atom_coordinates, best_response_value, interaction_manager):
        self.global_indices = global_indices
        self.type_array = type_array
        self.bonds = np.array(bonds)
        self.bond_orders = bond_orders
        self.renderer = v.add_renderer(BallAndStickRenderer, atom_coordinates[self.global_indices], type_array[global_indices], np.array(bonds_reduced))
        self.interaction_manager = interaction_manager

        """
        for i in self.global_indices:
            for j in self.global_indices:
                self.interaction_manager.interaction_matrix[i][j] = 9
                pass
        """

    def fix_carbon_bond_lengths(self, atom_coordinates):
        for bond in self.bonds:
            if bond[0] in self.global_indices and bond[1] in self.global_indices:
                bond_types = [self.type_array[bond[0]], self.type_array[bond[1]]]
                #print(bond, bond_types)
                if bond_types.count("C") == 2:
                    bond_length = 0.1 * self.interaction_manager.get_bond_length(*bond)
                    f, uf = self.bond_bisect(self.bonds, bond)
                    if np.linalg.norm(atom_coordinates[bond[1]] - atom_coordinates[bond[0]]) != bond_length:
                        ab = atom_coordinates[bond[1]] - atom_coordinates[bond[0]]
                        #print(("Distance from %s to %s"%(bond[0], bond[1]),"Before:", np.linalg.norm(ab)))
                        normalised = ab * (bond_length/np.linalg.norm(ab))
                        atom_coordinates[uf] +=(normalised-ab)
                        ab = atom_coordinates[bond[1]] - atom_coordinates[bond[0]]
                        #print(("Corrected", np.linalg.norm(ab)))
        return atom_coordinates

    def reduce(self, atom_coordinates):
        frag_distances = {}
        for frag_key in global_frag_distances.keys():
            if frag_key[0] - 1 in self.global_indices and frag_key[1] - 1 in self.global_indices:
                frag_distances[(frag_key[0] - 1, frag_key[1] - 1)] = global_frag_distances[frag_key]

        best_value = 10000
        best_atom_coordinates = atom_coordinates
        no_improve = 0

        while True:
            for x in range(0,100):
                atom_coordinates = self.rotate_bonds(best_atom_coordinates)
                print(atom_coordinates == best_atom_coordinates)
            value = self.get_noesy_value(frag_distances, atom_coordinates)
            if value < best_value-0.01:
                best_value = value
                best_atom_coordinates = np.copy(atom_coordinates)
                no_improve = 0
            else:
                no_improve += 1
            if no_improve > 1000:
                break

        return best_atom_coordinates

    def get_noesy_value(self, frag_distances, atom_coordinates):
        value = 0
        for frag_key in frag_distances.keys():

            vec = atom_coordinates[frag_key[0]] - atom_coordinates[frag_key[1]]

            if frag_key == (6,7):
                input("Distances for CHIRAL")
                input(np.linalg.norm(atom_coordinates[6] - atom_coordinates[14]))
                input(np.linalg.norm(atom_coordinates[7] - atom_coordinates[14]))


            input(("Frag Data", frag_key ,np.linalg.norm(10*vec), frag_distances[frag_key]))
            dist = np.linalg.norm(10*vec) - frag_distances[frag_key]
            value += dist * dist
        return value

    def update(self, best_atom_coordinates, best_response_value):
        #TODO REMOVE
        #best_atom_coordinates = self.reduce(best_atom_coordinates)

        for x in range(5):
            op = random.choice([self.rotate_bonds, self.translate])
            new_coordinates = op(best_atom_coordinates)
        new_response_value = calculate_response_value(new_coordinates, self.interaction_manager)

        if new_response_value < best_response_value:
            best_response_value = new_response_value
            best_atom_coordinates = new_coordinates
            self.renderer.update_positions(2*new_coordinates[self.global_indices])

        return best_atom_coordinates, best_response_value

    def rotate_bonds(self, atom_coordinates, bond=None):
        r_array = np.copy(atom_coordinates)
        if len(self.bonds) > 0:
            if bond is None:
                bond = random.choice(self.bonds)
            bond = self.bonds[0]
            """
            if 7 in bond or 6 in bond:
                input("6/7 rotation")
            """
            i1 = bond[0]
            i2 = bond[1]

            point_a = np.copy(atom_coordinates[i1])
            point_b = np.copy(atom_coordinates[i2])

            uf = random.choice(self.bond_bisect(self.bonds, bond, True))

            factor = 1000

            x = random.uniform(0, 2*factor)

            bond_vector = point_b - point_a
            random_vector = np.random.rand(3)
            rotation_axis = np.cross(bond_vector, random_vector)
            rotation = rotation_matrix(x*np.pi/factor, rotation_axis)[:3, :3]
            r_array -= point_a
            r_array[uf] = np.dot(r_array[uf], rotation.T)
            r_array += point_a
        else:
            i = random.choice(self.global_indices)
            j = random.randrange(3)
            r_array[i][j] += random.uniform(-0.1, 0.1)

        return r_array

    def translate(self, atom_coordinates):
        r_array = np.copy(atom_coordinates)
        r_array[self.global_indices] += random.uniform(-0.2, 0.2)
        return r_array

    def bond_bisect(self, bonds, bond, freeze_left=True):
        bonds = bonds.tolist()
        bond = bond.tolist()
        atoms = set()
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
                    if a == bond[0] and a not in atoms:
                        new_atoms.add(bond[1])
                    if a == bond[1] and a not in atoms:
                        new_atoms.add(bond[0])
            atoms.add(a)
            new_atoms.remove(a)
        frozen_atoms = atoms
        unfrozen_atoms = [x for x in list(set(sum(bonds, []))) if x not in frozen_atoms]
        return list(frozen_atoms), unfrozen_atoms


def calculate_response_value(atom_coordinates, interaction_manager):
        response = 0
        for i, v1 in enumerate(atom_coordinates):
            for j, v2 in enumerate(atom_coordinates):
                if j < i:
                    response += interaction_manager.interaction_response(i, j, v1*100, v2*100)
        return response

if __name__ == "__main__":
    clm = ChemLabMinimiser()
    clm.main()