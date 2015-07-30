import numpy as np
from chemlab.core import Atom, Molecule, System
from chemlab.graphics.renderers import BallAndStickRenderer
from chemlab.graphics import QtViewer
from chemlab.graphics.transformations import rotation_matrix
from chemlab.io.handlers import MolIO

import random
from main import get_interaction_manager,get_twod_signal_manager

class MolecularMinimizer:
    def __init__(self, molecules, interaction_manager):

        self.molecules = molecules
        for molecule in self.molecules:
            molecule.r_array *= 1.85
        self.v = QtViewer()
        self.renderers = {}
        for molecule in self.molecules:
            ar = self.v.add_renderer(BallAndStickRenderer, molecule.r_array, molecule.type_array, molecule.bonds)
            self.renderers[molecule] = ar
        self.best_response_value = 0
        self.best_atom_coordinates = None
        self.twod_signal_manager = get_twod_signal_manager()
        self.interaction_manager = interaction_manager

    def update(self):
        for molecule in self.molecules:
            if len(molecule.bonds)>=2:

                i = random.randrange(len(molecule.bonds))
                a = bond_bisect(molecule.bonds, molecule.bonds[i], True)

                r_old = np.copy(molecule.r_array)

                uf = a[1]
                point_a = np.copy(molecule.r_array[molecule.bonds[i][0]])
                point_b = np.copy(molecule.r_array[molecule.bonds[i][1]])
                M = rotation_matrix(np.pi/20, point_b-point_a)[:3, :3]

                molecule.r_array = np.array(molecule.export['vectors'])

                molecule.r_array -= point_a
                molecule.r_array[uf] = np.dot(molecule.r_array[uf], M.T)
                molecule.r_array += point_a
                self.set_vectors(molecule, molecule.r_array)
                new_response_value = calculate_response_value(molecule.r_array)
                if new_response_value < self.best_response_value:
                    self.best_response_value = new_response_value
                else:
                    molecule.r_array = r_old
                    self.set_vectors(r_old)
                self.renderers[molecule].update_positions(molecule.r_array)
                self.v.widget.repaint()

            else:
                None
    def set_vectors(self, molecule, arr):
        for i,vec in enumerate(arr):
            molecule.export['vectors'][i][:] = vec

    def run(self):
        self.v.schedule(self.update)
        self.v.run()

def split_groups(mol):
    groups = []
    unselected = [x for x in range(mol.n_atoms)]
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

            for bond_index, bond in enumerate(mol.bonds):
                if atoms[-1] == bond[0] and bond[1] not in atoms:
                    new_atoms.add(bond[1])
                if atoms[-1] == bond[1] and bond[0] not in atoms:
                    new_atoms.add(bond[0])
                if a in bond and bond.tolist() not in new_bonds:
                    new_bonds.append(bond.tolist())
                    new_bond_orders.append(mol.bond_orders[bond_index])
        atoms.sort()
        groups.append([atoms, new_bonds, new_bond_orders])
        unselected = [x for x in range(mol.n_atoms) if x not in sum([y[0] for y in groups],[])]
    molecules = []
    for group in groups:
        new_bonds_reduced = np.array([[group[0].index(x[0]), group[0].index(x[1])] for x in group[1]])
        numpy_vectors = [system.r_array[x] for x in group[0]]
        new_molecule = Molecule.from_arrays(r_array=np.array(numpy_vectors),
                                            type_array=system.type_array[a],
                                            bonds=new_bonds_reduced,
                                            bond_orders=group[2]
                                            )
        new_molecule.export['vectors'] = numpy_vectors
        molecules.append(new_molecule)
    return molecules

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

def calculate_response_value(atom_coordinates):
    response = 0
    for i, v1 in enumerate(atom_coordinates):
        for j, v2 in enumerate(atom_coordinates):
            if j < i:
                response += interaction_manager.interaction_response(i, j, v2, v1)
    return response


infile = open('../outfile.mol', 'rb')
molio = MolIO(infile)
system = molio.read('molecule')
twod_signal_manager = get_twod_signal_manager()
print(twod_signal_manager.get_interaction_matrix())
interaction_manager = get_interaction_manager(twod_signal_manager.get_interaction_matrix(), system.type_array.tolist(), np.zeros(shape=system.type_array.shape))

mm = MolecularMinimizer(split_groups(system), interaction_manager)
mm.run()
