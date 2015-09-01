import numpy as np


class CoordinateManager:
    def __init__(self, atom_coordinates):
        self.atom_coordinates = atom_coordinates
        self.old_atom_coordinates = np.copy(atom_coordinates)

    def start_iter(self):
        del self.old_atom_coordinates
        self.old_atom_coordinates = np.copy(self.atom_coordinates)

    def update(self, atom_coordinates):
        del self.atom_coordinates
        atom_coordinates -= atom_coordinates[-1]
        self.atom_coordinates = atom_coordinates

    def force_update(self, atom_coordinates):
        del self.atom_coordinates
        del self.old_atom_coordinates
        self.atom_coordinates = atom_coordinates
        self.old_atom_coordinates = np.copy(atom_coordinates)

    def revert(self):
        del self.atom_coordinates
        self.atom_coordinates = self.old_atom_coordinates

    def get_coordinates(self):
        return np.copy(self.atom_coordinates)

    def add_atom(self):
        self.atom_coordinates = np.insert(self.atom_coordinates, len(self.atom_coordinates), values=0.1, axis=0)
        del self.old_atom_coordinates
        self.old_atom_coordinates = np.copy(self.atom_coordinates)

