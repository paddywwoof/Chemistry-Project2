from chemlab.mviewer.qtmolecularviewer import QtMolecularViewer
from chemlab.mviewer.representations.ballandstick import BallAndStickRepresentation
from chemlab.core import System, Atom, Molecule
import numpy as np
from myrenderers.mybondrenderer import MyBondRenderer
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


class MyBallAndStickRepresentation(BallAndStickRepresentation):
    def __init__(self, interactions, interaction_colors, *args):
        super(MyBallAndStickRepresentation, self).__init__(*args)
        self.interaction_renderer = self.viewer.add_renderer(MyBondRenderer, interactions, interaction_colors, self.system.r_array, style='impostors', radius=0.005)

    def update_positions(self, r_array):
        super(MyBallAndStickRepresentation, self).update_positions(r_array)
        self.interaction_renderer.update_positions(r_array)

class MolecularGraphics:
    def __init__(self, interaction_manager):
        type_array = [x.atom_type for x in interaction_manager.atoms]
        bonds = interaction_manager.get_bonds_array()
        self.interaction_manager = interaction_manager
        self.coordinates = self.interaction_manager.get_coordinates()
        self.system = self.get_system()
        self.active = True
        print("Initialising Graphics Handler")
        self.viewer = QtMolecularViewer()
        self.update()
        if len(self.coordinates) != len(interaction_manager.atoms):
            raise Exception("Atom List does not correspond to coordinate array")

    def get_system(self):
        self.coordinates = self.interaction_manager.get_coordinates()
        atoms = []
        for index in range(len(self.interaction_manager.atoms)):
            atom = self.interaction_manager.atoms[index]
            type = atom.atom_type
            if type.upper() == "C":
                if atom.get_free_valency()!=0:
                    type = "Ne"
            pos = 3.5 * self.coordinates[index]
            atoms.append(Atom(type,pos))
        mol = Molecule(atoms, self.interaction_manager.get_bonds_array())
        system = System([mol])
        system.bonds = np.array(self.interaction_manager.get_bonds_array())
        return system

    def get_namespace(self):
        return self.viewer.namespace['__builtins__']

    def print_indices(self):
        indices = self.viewer.representation.selection_state['atoms'].indices
        if len(indices) > 0:
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

    def run(self, iterfunction=None):
        if self.active:
            if iterfunction:
                self.timer = self.viewer.schedule(iterfunction, 200)
            self.viewer.run()
        if not self.active:
            while True:
                iterfunction()

    def update_coordinates(self):
        self.viewer.representation.update_positions(3.5*self.interaction_manager.get_coordinates())

    def get_selected_atoms(self):
        selected_indices = self.viewer.representation.selection_state['atoms'].indices
        return [self.interaction_manager.atoms[x] for x in selected_indices]

    def update(self):
        self.system = self.get_system()
        self.viewer.clear()
        interactions, interaction_colors = self.get_through_space()
        self.viewer.system = self.system
        self.viewer.representation = MyBallAndStickRepresentation(interactions, interaction_colors,
                                                                  self.viewer, self.system)
