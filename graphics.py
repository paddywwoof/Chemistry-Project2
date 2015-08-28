
from chemlab.mviewer.qtmolecularviewer import QtMolecularViewer
from chemlab.mviewer.representations.ballandstick import BallAndStickRepresentation
from chemlab.core import System, Atom, Molecule
import numpy as np
from myrenderers.mybondrenderer import MyBondRenderer
from chemlab.graphics.renderers import BondRenderer
from signalmanager import InteractionValues



class Colors:
    red = (255, 0, 0, 255)
    orange = (255, 165, 0, 0)
    green = (0, 255, 0, 127)
    blue = (0, 0, 255, 255)
    purple = (128, 0, 255, 127)
    grey = (127, 127, 127, 255)
    pink = (255, 20, 147, 127)
    white = (255,255,255,255)


interaction_colors_dict = {InteractionValues.COSY: Colors.pink,
                           InteractionValues.HMBC: Colors.green,
                           InteractionValues.NOESY: Colors.purple}
bond_colors_dict = {2: Colors.orange, 3: Colors.blue, 4: Colors.grey, 7: Colors.red}


def get_MyBallAndStickRepresentation(interactions, interaction_colors):
    class MyBallAndStickRepresentation(BallAndStickRepresentation):
        def __init__(self, *args):
            super(MyBallAndStickRepresentation, self).__init__(*args)

            if len(interactions) != len(interaction_colors):
                raise Exception("Should be 1-1 correspondence between interaction and color")

            self.interaction_renderer = self.viewer.add_renderer(MyBondRenderer, interactions, interaction_colors, self.system.r_array, style='impostors', radius=0.005)

        def update_positions(self, r_array):
            super(MyBallAndStickRepresentation, self).update_positions(r_array)
            self.interaction_renderer.update_positions(r_array)
    return MyBallAndStickRepresentation


class MolecularGraphics:
    def __init__(self, coordinates, interaction_manager):
        type_array = [x.atom_type for x in interaction_manager.atoms]
        bonds = interaction_manager.get_bonds_array()
        self.interaction_manager = interaction_manager
        self.system = System.from_arrays(r_array=coordinates, type_array=type_array, mol_indices=[0])
        self.system.bonds = bonds
        self.active = True
        self.viewer = QtMolecularViewer()
        self.coordinates = coordinates
        self.system = self.get_system()

        self.viewer.add_representation(get_MyBallAndStickRepresentation(*self.get_through_space()), self.system)

        try:
            self.update(coordinates)
            if len(coordinates) != len(interaction_manager.atoms):
                raise Exception("Atom List does not correspond to coordinate array")
        except Exception as e:
            self.disable(e)

    def get_system(self):
        atoms = []
        for index in range(len(self.interaction_manager.atoms)):
            type = self.interaction_manager.atoms[index].atom_type
            pos = self.coordinates[index]
            atoms.append(Atom(type,pos))
        mol = Molecule(atoms, self.interaction_manager.get_bonds_array())
        system = System([mol])
        return system

    def add_bondsetter(self, bond_setter):
        self.viewer.namespace['__builtins__'].addbond = bond_setter

    def add_atomsetter(self, atom_setter):
        self.viewer.namespace['__builtins__'].addatom = atom_setter

    def print_indices(self):
        indices = self.viewer.representation.selection_state['atoms'].indices
        if len(indices)>0:
            print(indices)

    def disable(self, e):
        input("Exception occured in graphics handler. Graphics will henceforth be disabled:")
        self.active = False
        raise e


    def get_through_space(self):
        interactions = []
        interaction_colors = []

        for interaction_number in InteractionValues.THROUGH_SPACE:
            new_interaction_list = self.interaction_manager.get_all_interaction_atoms(interaction_number)
            new_interaction_indices = [[y[0].index_value, y[1].index_value] for y in new_interaction_list]

            interactions += new_interaction_indices
            interaction_colors += [interaction_colors_dict[interaction_number]]*len(new_interaction_indices)
        return np.array(interactions), np.array(interaction_colors)


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
        self.coordinates = coordinates
        self.system = self.get_system()
        self.system.bonds = np.array(self.interaction_manager.get_bonds_array())
        self.viewer.clear()
        self.viewer.add_representation(get_MyBallAndStickRepresentation(*self.get_through_space()), self.system)
