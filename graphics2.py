
from chemlab.mviewer.qtmolecularviewer import QtMolecularViewer
from chemlab.mviewer.representations.ballandstick import BallAndStickRepresentation
from chemlab.core import System, Atom, Molecule
import numpy as np

from chemlab.graphics.renderers import BondRenderer


class MyBallAndStickRepresentation(BallAndStickRepresentation):
    def on_atom_selection_changed(self):
        sel = self.selection_state['atoms'].indices
        cols = np.array(self.atom_colors.array)
        cols[sel, -1] = 50
        self.atom_renderer.update_colors(cols)
        self.viewer.update()

    def update_bond_renderer(self):
        self.viewer.remove_renderer(self.bond_renderer)
        self.bond_renderer = self.viewer.add_renderer(BondRenderer,
                                                      self.system.bonds, self.system.r_array,
                                                      self.system.type_array, style='impostors')
        self.viewer.update()

class MolecularGraphics:
    def __init__(self, coordinates, interaction_manager):
        type_array = [x.atom_type for x in interaction_manager.atoms]
        bonds = interaction_manager.get_bonds_array()
        self.interaction_manager = interaction_manager
        self.system = System.from_arrays(r_array=coordinates, type_array=type_array, mol_indices=[0])
        self.system.bonds = bonds
        self.active = True
        self.viewer = QtMolecularViewer()

        atoms = []
        for index in range(len(coordinates)):
            type = interaction_manager.atoms[index].atom_type
            pos = coordinates[index]
            atoms.append(Atom(type,pos))
        mol = Molecule(atoms, interaction_manager.get_bonds_array())

        self.system = System([mol])
        self.viewer.add_representation(MyBallAndStickRepresentation, self.system)



        try:
            self.update(coordinates)
            if len(coordinates) != len(interaction_manager.atoms):
                raise Exception("Atom List does not correspond to coordinate array")
        except Exception as e:
            self.disable(e)

    def print_indices(self):
        indices = self.viewer.representation.selection_state['atoms'].indices
        if len(indices)>0:
            print(indices)

    def disable(self, e):
        input("Exception occured in graphics handler. Graphics will henceforth be disabled:")
        self.active = False
        raise e

    def run(self, iterfunction):
        if self.active:
            try:
                self.viewer.schedule(iterfunction)
                self.viewer.run()
            except Exception as e:
                self.disable(e)
        if not self.active:
            while True:
                iterfunction()

    def update_coordinates(self, coordinates):
        self.viewer.representation.update_positions(3.5*coordinates)

    def get_selected_atoms(self):
        selected_indices = self.viewer.representation.selection_state['atoms'].indices
        return [self.interaction_manager.atoms[x] for x in selected_indices]

    def update(self, coordinates):
        self.system.r_array = 3.5*coordinates
        self.system.bonds = np.array(self.interaction_manager.get_bonds_array())
        self.viewer.clear()
        self.viewer.add_representation(MyBallAndStickRepresentation, self.system)

