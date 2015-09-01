import numpy as np


class SystemState:
    def __init__(self, message, atoms, bonds, coordinates, interaction_matrix):
        self.message = message
        self.atom_data = [atom.export() for atom in atoms]
        self.bond_data = [bond.export() for bond in bonds]
        self.coordinates = np.copy(coordinates)
        self.interaction_matrix = np.copy(interaction_matrix)
        self.Atom = atoms[0].__class__
        self.Bond = bonds[0].__class__

    def get_data(self):
        atoms = []
        for x in self.atom_data:
            new_atom = self.Atom(x[0], x[1], x[2])
            atoms.append(new_atom)

        bonds = []
        for x in self.bond_data:
            atom1 = atoms[x[0][0]]
            atom2 = atoms[x[0][1]]
            new_bond = self.Bond(atom1, atom2, x[1], x[2], x[3])
            bonds.append(new_bond)
        return atoms, bonds
