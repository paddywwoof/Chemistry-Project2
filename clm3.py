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
        im = self.interaction_manager.interaction_matrix
        im[19][6] = 0
        im[19][7] = 0

        """
        for x in [9,10,11,12,13,17,18,19,20,21]:
            for y in range(self.interaction_manager.number_atoms):
                if im[x][y] == 3:
                    im[x][y] = 0
                    im[y][x] = 0

        """

        self.interaction_matrix = self.interaction_manager.interaction_matrix
        self.type_array = self.interaction_manager.type_array
        self.shift_data = self.interaction_manager.shift_data
        self.number_atoms = self.interaction_manager.number_atoms
        self.global_bond_indices = self.interaction_manager.bonds
        self.bond_orders = self.interaction_manager.bond_orders
        self.viewer = QtViewer()
        self.best_response_value = 1000000.0
        ac = np.random.rand(self.number_atoms, 3)*0.001
        #ac = FileManager().read_numpy_from_xyz('manual2.xyz') * 0.1
        self.iteration_number = 0
        self.coordinate_manager = CoordinateManager(self, ac, self.interaction_manager)
        self.fragments = self.generate_fragments_list()
        print(self.coordinate_manager.calculate_response_value(True))
        self.interaction_manager.print_matrix()
        input("")

        self.delay = 0
        """
        for x in [0,1,2,3,4,5]:
            self.interaction_manager.interaction_matrix[self.interaction_manager.interaction_matrix==x] = 9
        """
        self.nullfrag = Fragment([], [], [], [], self.viewer, self.coordinate_manager, self.interaction_manager, [])




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
            local_bond_indices = np.array([[global_indices.index(x[0]), global_indices.index(x[1])] for x in global_bond_indices])
            bond_orders = fragment_piece['bond orders']
            fragment = Fragment(global_indices,
                                global_bond_indices,
                                local_bond_indices,
                                bond_orders,
                                self.viewer,
                                self.coordinate_manager,
                                self.interaction_manager,
                                self.free_valencies
                                )
            fragments.append(fragment)
        return fragments

    def calculate_valencies(self):
        valencies = {"H":1, "C":4, "N":4}
        self.free_valencies = [valencies[x] for x in self.interaction_manager.type_array]
        for bond in self.global_bond_indices:
            self.free_valencies[bond[0]]-=1
            self.free_valencies[bond[1]]-=1
        for i, shift in enumerate(self.shift_data):
            if 100 < shift < 160:
                self.free_valencies[i]-=1
        self.fix_types()

    def fix_types(self):
        for index, type in enumerate(self.type_array):
            if self.free_valencies[index] == 0 and self.type_array[index] == "C":
                self.type_array[index]="N"

    def fix_renderers(self):
        self.viewer.clear()
        for fragment in self.fragments:
            fragment.set_renderer()



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
        #self.interaction_manager.interaction_matrix[self.interaction_manager.interaction_matrix==1] = 0
        #self.interaction_manager.interaction_matrix[self.interaction_manager.interaction_matrix==2] = 0

        print(self.coordinate_manager.calculate_response_value(True))
        self.viewer.schedule(self.iteration2)
        self.viewer.run()

    def iteration2(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        hmbcs = []
        for x1 in range(len(atom_coordinates)):
            for x2 in range(len(atom_coordinates)):
                if self.interaction_matrix[x1][x2]==3:
                    hmbcs.append([x1,x2])
        for frag1 in self.fragments:
            for frag2 in self.fragments:
                for hmbc in hmbcs:
                    a = hmbc[0] in frag1.global_indices and hmbc[1] in frag2.global_bond_indices
                    b = hmbc[1] in frag1.global_indices and hmbc[0] in frag2.global_bond_indices
                    if not (a and b):
                        continue
                    hmbc.sort()
                    if self.type_array[hmbc[0]] == "H":
                        bond = [x for x in self.global_bond_indices if hmbc[0] in x][0]
                        carbon_index = max(bond)
                        if self.free_valencies[carbon_index]!=0:
                            continue
                        else:
                            adjacent_carbons = [[y for y in x if y!=carbon_index and self.free_valencies[y]!=0][0] for x in self.global_bond_indices if carbon_index in x]
                            if len(adjacent_carbons) == 1:
                                self.merge(frag1,frag2,)
                                
                                
                                self.global_bond_indices.append(adjacent_carbons[0], hmbc[1])

    def iteration(self):
        self.iteration_number+=1
        if self.iteration_number < 2500:
            self.optimise_substructures()
        else:
            self.form_new_bonds()

    def optimise_substructures(self):
        for fragment in self.fragments:
            fragment.start_iter()
            fragment.rotate_bonds(sequence=1)
            fragment.rotate_bonds(sequence=1)
            fragment.verify()
            fragment.start_iter()
            fragment.translate()
            fragment.verify()

            fragment.update_graphics()

    def form_new_bonds(self):


        fragment1 = random.choice(self.fragments)

        shortest_bond = None
        closest_fragment = None
        shortest_distance = 1000

        atom_coordinates = self.coordinate_manager.get_coordinates()

        for fragment2 in self.fragments:
            if fragment1 != fragment2 or len(self.fragments)==1:
                for i1 in fragment1.global_indices:
                    for i2 in fragment2.global_indices:
                        frag1_atom = atom_coordinates[i1]
                        frag2_atom = atom_coordinates[i2]
                        if np.linalg.norm(frag2_atom-frag1_atom) < shortest_distance:
                            if self.free_valencies[i1]>0 and self.free_valencies[i2]>0 and self.checkHMBC(i1, i2, fragment1, fragment2):
                                if not(fragment1.has_hydrogens[i1] and fragment1.has_hydrogens[i2]):
                                    shortest_bond = [i1, i2]
                                    shortest_distance = np.linalg.norm(frag2_atom-frag1_atom)
                                    closest_fragment = fragment2

        if closest_fragment == None:
            return
        else:
            self.iteration_number=0
            self.coordinate_manager.reset()


        self.interaction_manager.interaction_matrix[i1][i2]=4
        self.interaction_manager.interaction_matrix[i2][i1]=4
        if len(self.fragments)>1:
            self.merge(fragment1, closest_fragment, shortest_bond)
        else:


            self.merge(fragment1, self.nullfrag, shortest_bond)


        self.coordinate_manager.reset()
        for fragment in self.fragments:
            fragment.project_bond_lengths()


    def merge(self, fragment1, fragment2, shortest_bond):
        """
        Creates a list of fragment objects
        """

        global_indices = fragment1.global_indices + fragment2.global_indices
        self.global_bond_indices.append(shortest_bond)
        self.calculate_valencies()
        global_bond_indices = fragment1.global_bond_indices + fragment2.global_bond_indices + [shortest_bond]
        local_bond_indices = np.array([[global_indices.index(x[0]), global_indices.index(x[1])] for x in global_bond_indices])
        bond_orders = fragment1.bond_orders + fragment2.bond_orders


        self.fix_types()
        self.fix_renderers()

        fragment = Fragment(global_indices,
                            global_bond_indices,
                            local_bond_indices,
                            bond_orders,
                            self.viewer,
                            self.coordinate_manager,
                            self.interaction_manager,
                            self.free_valencies
                            )
        self.viewer.clear()
        if fragment1 in self.fragments:
            self.fragments.remove(fragment1)
        if fragment2 in self.fragments:
            self.fragments.remove(fragment2)
        self.fragments.append(fragment)
        self.calculate_valencies()
        for fragment in self.fragments:
            fragment.set_renderer()

        return fragment

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
    def __init__(self, global_indices, global_bond_indices, local_bond_indices, bond_orders, viewer, coordinate_manager, interaction_manager, free_valencies):
        self.global_indices = global_indices
        self.global_bond_indices = global_bond_indices
        self.local_bond_indices = local_bond_indices
        self.bond_orders = bond_orders
        self.type_array = np.array(interaction_manager.type_array)
        self.coordinate_manager = coordinate_manager
        self.viewer = viewer
        self.coordinate_manager = coordinate_manager
        self.interaction_manager = interaction_manager
        self.project_bond_lengths()
        self.number_atoms = len(self.global_indices)
        self.best_response_value = 1000
        self.verify_index = 0
        self.set_renderer()
        self.free_valencies = free_valencies

        self.has_hydrogens = [False for x in range(len(self.coordinate_manager.atom_coordinates))]
        for x in self.global_bond_indices:
            if self.type_array[x[0]]=="H":
                self.has_hydrogens[x[1]]=True
            if self.type_array[x[1]]=="H":
                self.has_hydrogens[x[0]]=True

    def add_bond(self, bond):
        if bond[0] not in self.global_indices and bond[1] not in self.global_indices:
            raise Exception
        else:
            if bond[1] in self.global_indices:
                bond.reverse()

        self.global_indices.append(bond[1])
        local_bond_index_1 = self.global_indices.index(bond[0])
        local_bond_index_2 = self.global_indices.index(bond[1])

        self.global_bond_indices.append(bond)
        self.number_atoms = len(self.global_indices)
        self.local_bond_indices = np.append(self.local_bond_indices,[local_bond_index_1,local_bond_index_2])


    def project_bond_lengths(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for bond in self.global_bond_indices:
            if self.interaction_manager.interaction_matrix[bond[0]][bond[1]]==9:
                continue
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
        unfrozen_atoms = [x for x in list(set(sum(bonds, []))) if x not in frozen_atoms]#


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

        xvec = (np.random.rand(3)-0.5)
        xvec /= np.linalg.norm(xvec)

        vec = random.uniform(0,0.6)*xvec
        new_atom_coordinates[self.global_indices] += vec

        self.coordinate_manager.update(new_atom_coordinates)

    def rotate(self):
        new_atom_coordinates = self.coordinate_manager.get_coordinates()
        atom_index = random.choice(self.global_indices)
        atom = new_atom_coordinates[atom_index]
        new_atom_coordinates -= atom
        angle = random.uniform(0, np.pi)
        rotation_axis = np.random.rand(3)
        rotation = rotation_matrix(angle, rotation_axis)[:3, :3]
        new_atom_coordinates[self.global_indices] = np.dot(new_atom_coordinates[self.global_indices], rotation.T)
        new_atom_coordinates += atom
        self.coordinate_manager.update(new_atom_coordinates)
    def update_graphics(self):
        x = np.random.normal(loc=3.0, scale=3.0)

        self.graphics.update_positions(2*self.coordinate_manager.atom_coordinates[self.global_indices])




class CoordinateManager:
    def __init__(self, parent, atom_coordinates, interaction_manager):
        self.parent = parent
        self.atom_coordinates = atom_coordinates
        self.old_atom_coordinates = np.copy(atom_coordinates)

        self.interaction_manager = interaction_manager
        self.best_response_value = self.calculate_response_value()
        self.file_manager = FileManager()
        self.ival_matrix = np.zeros_like(self.interaction_manager.interaction_matrix, dtype=float)
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
        self.temperature /= 1.001



        response = self.calculate_response_value()
        deltaE = response - self.best_response_value
        if deltaE <= 0:
            self.best_response_value = response
            print("New Best Response: %s at iteration %s"%(response,self.parent.iteration_number))
            self.calculate_response_value(True)
            self.file_manager.write_numpy_to_mol("tempfile.mol", self.interaction_manager, self.atom_coordinates)
            return True
        else:
            self.revert()
            return False

    def revert(self):
        self.atom_coordinates = self.old_atom_coordinates

    def calculate_response_value(self, debug=False):
        response = 0
        if debug:
            print("_*_"*3)
        for i, v1 in enumerate(self.atom_coordinates):
            for j, v2 in enumerate(self.atom_coordinates):
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, v1, v2, debug)
        return response

    def get_coordinates(self):
        return np.copy(self.atom_coordinates)

def end():
    input("Press Any Key to Quit")
    raise SystemExit

def main():
    clm = ChemLabMinimiser()
    file_manager = FileManager()
    clm.main()
    end()

if __name__ == "__main__":
    main()
