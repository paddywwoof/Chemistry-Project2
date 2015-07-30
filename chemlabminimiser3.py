from main_run import get_interaction_manager,get_twod_signal_manager
from filemanager import FileManager
from chemlab.graphics import QtViewer
import numpy as np
import random
from chemlab.graphics.transformations import rotation_matrix
from chemlab.graphics.renderers import BallAndStickRenderer
from math import exp
import os

class ChemLabMinimiser:
    def __init__(self):
        self.interaction_manager = get_interaction_manager(*get_twod_signal_manager().get_interaction_data())
        self.interaction_manager.interaction_matrix[self.interaction_manager.interaction_matrix==1] = 0
        self.interaction_manager.interaction_matrix[self.interaction_manager.interaction_matrix==3] = 0


        self.interaction_matrix = self.interaction_manager.interaction_matrix



        self.type_array = self.interaction_manager.type_array
        self.shift_data = self.interaction_manager.shift_data
        self.number_atoms = self.interaction_manager.number_atoms
        self.global_bond_indices = self.interaction_manager.bonds
        self.bond_orders = self.interaction_manager.bond_orders
        self.viewer = QtViewer()
        self.best_response_value = 1000000.0
        ac = np.random.rand(self.number_atoms, 3)*3
        #ac = FileManager().read_numpy_from_xyz('manual2.xyz') * 0.1
        self.coordinate_manager = CoordinateManager(ac, self.interaction_manager)
        self.fragments = self.generate_fragments_list()
        print(self.coordinate_manager.calculate_response_value())


        self.delay = 0
        """
        for x in [0,1,2,3,4,5]:
            self.interaction_manager.interaction_matrix[self.interaction_manager.interaction_matrix==x] = 9
        """


    def generate_fragments_list(self):
        """
        Creates a list of fragment objects
        """
        fragments = []
        fragment_data = self.generate_fragment_groups()
        for fragment_piece in fragment_data:
            global_indices = fragment_piece['global indices']
            global_bond_indices = fragment_piece['global bond indices']
            local_bond_indices = np.array([[global_indices.index(x[0]), global_indices.index(x[1])] for x in global_bond_indices])
            bond_orders = fragment_piece['bond orders']
            fragment = Fragment(global_indices,
                                global_bond_indices,
                                local_bond_indices,
                                bond_orders,
                                self.viewer,
                                self.coordinate_manager,
                                self.interaction_manager
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
                for bond_index, bond in enumerate(self.global_bond_indices):
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
        for fragment in self.fragments:
            fragment.update_graphics()
        self.viewer.schedule(self.iteration)
        self.viewer.run()

    def iteration(self):
        if self.delay > 0:
            self.delay -=1
        if self.delay == 0:
            self.optimise_substructures()

    def optimise_substructures(self):
        for fragment in self.fragments:
            fragment.rotate_bonds(sequence=1)
            fragment.translate()
            fragment.rotate()
            fragment.update_graphics()


class Fragment:
    def __init__(self, global_indices, global_bond_indices, local_bond_indices, bond_orders, viewer, coordinate_manager, interaction_manager):
        self.global_indices = global_indices
        self.global_bond_indices = global_bond_indices
        self.local_bond_indices = local_bond_indices
        self.bond_orders = bond_orders
        type_array = np.array(interaction_manager.type_array)
        self.renderer = viewer.add_renderer(BallAndStickRenderer, coordinate_manager.atom_coordinates[global_indices], type_array[global_indices], np.array(local_bond_indices))
        self.coordinate_manager = coordinate_manager
        self.interaction_manager = interaction_manager
        self.project_bond_lengths()

        self.best_response_value = 1000
        self.verify_index = 0

    def project_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.global_bond_indices:
            bond_length = self.interaction_manager.get_bond_length(*bond)
            print(bond_length)
            frozen, unfrozen = self.bisect_on_bond(bond)
            bond_vector = atom_coordinates[bond[1]] - atom_coordinates[bond[0]]
            corrected_bond_vector = bond_vector * (bond_length/np.linalg.norm(bond_vector))
            atom_coordinates[unfrozen] += (corrected_bond_vector-bond_vector)
        self.coordinate_manager.update(atom_coordinates)

    def local_response_value(self):
        atom_coordinates = self.coordinate_manager.atom_coordinates
        response = 0
        for i in self.global_indices:
            for j in self.global_indices:
                v1 = atom_coordinates[i]
                v2 = atom_coordinates[j]
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, v1, v2)
        return response

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
        unfrozen_atoms = [x for x in list(set(sum(bonds, []))) if x not in frozen_atoms]
        return list(frozen_atoms), unfrozen_atoms

    def rotate_bonds(self, bond_index=None, angle=None, sequence=1):
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
            angle = (np.pi/100) * random.uniform(0, 200)

            if random.randrange(0,2)==0:
                rotation_axis = np.random.rand(3)
            else:
                rotation_axis = point_b - point_a

            rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
            new_atom_coordinates -= point_b
            new_atom_coordinates[uf] = np.dot(new_atom_coordinates[uf], rotation.T)
            new_atom_coordinates += point_b
            self.coordinate_manager.update(new_atom_coordinates)
            self.verify()
        else:
            self.translate()

    def verify(self):
        self.coordinate_manager.verify_global()


    def translate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()

        new_atom_coordinates[self.global_indices] += 0.1*(np.random.rand(3)-0.5)
        self.coordinate_manager.update(new_atom_coordinates)
        valid = self.verify()

    def rotate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        atom_index = random.choice(self.global_indices)
        atom = new_atom_coordinates[atom_index]
        new_atom_coordinates -= atom
        angle = random.uniform(0, 2*np.pi)
        rotation_axis = np.random.rand(3)
        rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
        new_atom_coordinates[self.global_indices] = np.dot(new_atom_coordinates[self.global_indices], rotation.T)
        new_atom_coordinates += atom
        self.coordinate_manager.update(new_atom_coordinates)
        self.verify()

    def update_graphics(self):
        x = np.random.normal(loc=3.0, scale=3.0)
        self.renderer.update_positions(2*self.coordinate_manager.atom_coordinates[self.global_indices])


class CoordinateManager:
    def __init__(self, atom_coordinates, interaction_manager):
        self.atom_coordinates = atom_coordinates
        self.old_atom_coordinates = np.copy(atom_coordinates)
        self.best_response_value = 1000000
        self.interaction_manager = interaction_manager
        self.file_manager = FileManager()
        self.ival_matrix = np.zeros_like(self.interaction_manager.interaction_matrix)

    def reset(self):
        self.best_response_value = 1000000

    def update(self, atom_coordinates):
        for x in [0, 1, 2]:
            atom_coordinates[:, x] -= min(atom_coordinates[:, x])

        self.old_atom_coordinates = np.copy(self.atom_coordinates)
        self.atom_coordinates = atom_coordinates

    def verify_global(self):
        response = self.calculate_response_value()
        if response < self.best_response_value:
            self.best_response_value = response
            print("New Best Response: %s"%response)
            self.calculate_response_value(True)
            self.file_manager.write_numpy_to_mol("tempfile.mol", self.interaction_manager, 1.2*self.atom_coordinates)
            return True
        else:
            self.revert()
            return False

    def revert(self):
        self.atom_coordinates = self.old_atom_coordinates

    def calculate_response_value(self, debug=False):
        new_ival_matrix = np.zeros_like(self.interaction_manager.interaction_matrix)
        response = 0
        for i, v1 in enumerate(self.atom_coordinates):
            for j, v2 in enumerate(self.atom_coordinates):
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, v1, v2, debug)
        return response

    def get_coordinates(self):
        return np.copy(self.atom_coordinates)


def main():
    clm = ChemLabMinimiser()
    file_manager = FileManager()
    clm.main()
    raise SystemExit

if __name__ == "__main__":
    main()
