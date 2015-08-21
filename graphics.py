from chemlab.graphics import QtViewer
#from chemlab.graphics.renderers import BallAndStickRenderer
from myrenderers.myballandstickrenderer import MyBallAndStickRenderer
import numpy as np
from chemlab.graphics.postprocessing import OutlineEffect

scale = 4




red = (255, 0, 0, 255)
orange = (255, 165, 0, 0)
green = (0,255,0,127)
blue = (0, 0, 255, 255)
purple = (128,0,255,127)

grey = (127, 127, 127, 255)
pink = (255,20,147,127)

Purple = (230, 101, 243)

interaction_colors_dict = {1: pink, 3: green, 6: purple}
bond_colors_dict = {2: orange, 3: blue, 4: grey, 7: red}


class MolecularGraphics:
    def __init__(self, coordinates, interaction_manager):
        self.viewer = QtViewer()
        self.active = True
        self.coordinates = None
        self.atoms = None
        self.bonds = None
        self.interaction_manager = interaction_manager
        try:
            self.update(coordinates, interaction_manager)
            if len(coordinates) != len(self.atoms):
                raise Exception("Atom List does not correspond to coordinate array")

        except Exception as e:
            self.disable(e)

    def update_coordinates(self, coordinates):
        if self.active:
            try:

                self.renderer.update_positions(scale*coordinates)
                self.viewer.widget.update()
            except Exception as e:
                self.disable(e)

    def disable(self, e):
        raise e
        input("Exception occured in graphics handler. Graphics will henceforth be disabled:")
        self.active = False

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

    def update(self, coordinates, interaction_manager):
        print("UPDATING GRAPHICS")
        self.atoms = interaction_manager.atoms
        self.bonds = interaction_manager.bonds
        if self.active:
            try:
                bond_colors = np.array([bond_colors_dict[bond.inferred_by] for bond in self.bonds])
                interactions = []
                interaction_colors = []
                for x in [1, 3, 6]:
                    new_interaction_list = interaction_manager.get_all_interaction_atoms(x)
                    new_interaction_indices = [[y[0].index_value, y[1].index_value] for y in new_interaction_list]

                    interactions += new_interaction_indices
                    interaction_colors += [interaction_colors_dict[x]]*len(new_interaction_indices)

                self.viewer.clear()
                self.renderer = self.viewer.add_renderer(MyBallAndStickRenderer,
                                                         scale * coordinates,
                                                         np.array([atom.atom_type for atom in self.atoms]),

                                                         np.array([bond.get_indices() for bond in self.bonds]),
                                                         bond_colors,

                                                         np.array(interactions),
                                                         np.array(interaction_colors)
                                                         )
            except Exception as e:
                self.disable(e)
