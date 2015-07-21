import numpy as np
from chemlab.graphics.renderers import BallAndStickRenderer
from chemlab.graphics import QtViewer
from chemlab.graphics.transformations import rotation_matrix
from chemlab.io.handlers import MolIO
import random

from main import get_interaction_manager,get_twod_signal_manager

def split_groups(system):
    #mol : mol.n_atoms, mol.bonds, mol.bond_orders
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

def bond_bisect(bonds, bond, freeze_left):
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


class Fragment:
    def __init__(self, global_indices, type_array, bonds, bonds_reduced, bond_orders, v, atom_coordinates, best_response_value):
        self.global_indices = global_indices
        self.type_array = type_array
        self.bonds = np.array(bonds)
        self.bond_orders = bond_orders
        self.atom_coordinates = atom_coordinates
        self.renderer = v.add_renderer(BallAndStickRenderer, self.atom_coordinates[self.global_indices], type_array[global_indices], np.array(bonds_reduced))
        self.best_response_value = best_response_value

        for i in self.global_indices:
            for j in self.global_indices:
                interaction_manager.interaction_matrix[i][j] = 9
                pass

    def check_new(self, r_array):
        new_response_value = calculate_response_value(r_array)
        if new_response_value < self.best_response_value[0]:
            self.best_response_value[0] = new_response_value
            self.best_atom_coordinates = r_array
            self.renderer.update_positions(2*r_array[self.global_indices])
            for i, vec in enumerate(r_array):
                self.atom_coordinates[i][:] = r_array[i]
            print(self.best_response_value[0])
            return True
        else:
            return False

    def update(self):

        updated = True
        while updated:
            r_array = self.rotate_bonds()
            updated = self.check_new(r_array)
        updated = True
        """
        while updated:
            r_array = self.translate()
            updated = self.check_new(r_array)

        updated = True
        while updated:
            r_array = self.atom_coordinates * 1.01
            updated = self.check_new(r_array)
        """


    def translate(self):
        r_array = np.copy(self.atom_coordinates)
        r_array[self.global_indices] += 0.05*np.random.rand(3)
        return r_array

    def rotate_bonds(self):
        r_array = np.copy(self.atom_coordinates)
        if len(self.bonds) > 0:
            i = random.randrange(len(self.bonds))

            i1 = self.bonds[i][0]
            i2 = self.bonds[i][1]
            point_a = np.copy(self.atom_coordinates[i1])
            point_b = np.copy(self.atom_coordinates[i2])

            uf = bond_bisect(self.bonds, self.bonds[i], True)[0]
            x = random.uniform(0, 40)

            bond_vector = point_b - point_a
            random_vector = np.random.rand(3)
            rotation_axis = np.cross(bond_vector, random_vector)


            M = rotation_matrix(x*np.pi/20, rotation_axis)[:3, :3]


            r_array -= point_a
            r_array[uf] = np.dot(r_array[uf], M.T)
            r_array += point_a
        else:

            i = random.choice(self.global_indices)
            j = random.randrange(3)
            """
            input(self.bonds)
            i, j = tuple(r_array.shape)
            i, j = (random.randrange(i), random.randrange(j))
            """
            r_array[i][j] += random.uniform(-3, 3)
        return r_array


def calculate_response_value(atom_coordinates):
    response = 0
    for i, v1 in enumerate(atom_coordinates):
        for j, v2 in enumerate(atom_coordinates):
            if j < i:
                response += interaction_manager.interaction_response(i, j, v1*100/scale, v2*100/scale)
    return response

def update():
    for x in [0,1,2]:
        atom_coordinates[:,x]-=min(atom_coordinates[:,x])
    for fragment in fragments:
        fragment.update()


    interaction_manager.bonds = system.bonds+1
    file_manager.write_numpy_to_mol("../resources/tempfile.mol", interaction_manager, atom_coordinates*10)


if __name__ == "__main__":
    from filemanager import FileManager
    file_manager = FileManager()
    scale = 1
    infile = open('../outfile.mol', 'rb')
    molio = MolIO(infile)
    system = molio.read('molecule')
    atom_coordinates = system.r_array * 1.818

    interaction_manager = get_interaction_manager(get_twod_signal_manager().get_interaction_matrix(), system.type_array.tolist(), np.zeros(shape=system.type_array.shape))
    #interaction_manager.interaction_matrix[interaction_manager.interaction_matrix!=3] = 9

    interaction_manager.print_matrix()
    interaction_manager.plot_all_interactions()
    groups = split_groups(system)
    v = QtViewer()
    fragments = []
    brv = np.array([1000000.0])



    for bond in system.bonds:
        bond_types = [system.type_array[bond[0]], system.type_array[bond[1]]]

        if bond_types == ["C", "C"]:
            print("CC Bond: ", np.linalg.norm(atom_coordinates[bond[0]] - atom_coordinates[bond[1]]))

        if "H" in bond_types and "C" in bond_types:
            if bond_types[0] == "H":
                h,c = bond
            else:
                c,h = bond
            point_a = np.copy(atom_coordinates[c])
            point_b = np.copy(atom_coordinates[h])

            vec = (point_b-point_a)
            vec /= np.linalg.norm(vec)
            vec *= scale * 0.11
            atom_coordinates[h][:] = point_a + vec



    for group in groups:
        new_bonds_reduced = np.array([[group[0].index(x[0]), group[0].index(x[1])] for x in group[1]])
        fragment = Fragment(group[0], system.type_array, np.array(group[1]),new_bonds_reduced, group[2], v, atom_coordinates, brv)
        fragments.append(fragment)
    v.schedule(update, 1)
    v.run()
raise SystemExit
