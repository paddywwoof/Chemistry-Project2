from chemlab.graphics import QtViewer
#from chemlab.graphics.renderers import BallAndStickRenderer
from myrenderers.myballandstickrenderer import MyBallAndStickRenderer
import numpy as np
#from chemlab.graphics.postprocessing import OutlineEffect
from chemlab.graphics.colors import default_atom_map

scale = 4

red = (255, 0, 0, 255)
orange = (255, 165, 0, 0)
green = (0, 255, 0, 127)
blue = (0, 0, 255, 255)
purple = (128, 0, 255, 127)

grey = (127, 127, 127, 255)
pink = (255, 20, 147, 127)
white = (255,255,255,255)


interaction_colors_dict = {1: pink, 3: green, 6: purple}
bond_colors_dict = {2: orange, 3: blue, 4: grey, 7: red}



atom_colors_dict = {"C": orange}


class MolecularGraphics:
    def __init__(self, coordinates, interaction_manager):
        self.viewer = QtViewer()
        self.selected_atom = 0
        self.active = True
        self.coordinates = None
        self.atoms = None
        self.bonds = None
        self.interaction_manager = interaction_manager

        self.viewer.widget.on_mouse.append(self.select)

        try:
            self.update(coordinates, interaction_manager)
            if len(coordinates) != len(self.atoms):
                raise Exception("Atom List does not correspond to coordinate array")

        except Exception as e:
            self.disable(e)

    def select(self, mouse_vector):
        closest_index = None
        closest_dist = 100
        for index, atom in enumerate(self.coordinates):
            dist = np.linalg.norm(mouse_vector[:2]-atom[:2])
            if dist < closest_dist:
                closest_index = index
                closest_dist = dist

        mouse_vector[2] = 0
        self.selected_atom = closest_index
        self.coordinates[0] = mouse_vector
        print(mouse_vector)
        self.update(self.coordinates, self.interaction_manager)



    """
    def scroll_atom(self):
        self.selected_atom += 1
        self.selected_atom %= len(self.interaction_manager.atoms)
        #print(self.selected_atom)
        self.update(self.coordinates, self.interaction_manager)
    """

    def update_coordinates(self, coordinates):
        self.coordinates = coordinates
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
        self.coordinates = coordinates
        self.interaction_manager = interaction_manager
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



                type_array = [atom.atom_type for atom in interaction_manager.atoms]
                atom_colors = []
                for index, atom in enumerate(interaction_manager.atoms):
                    if atom.get_free_valency() != 0 and atom.atom_type in atom_colors_dict and self.selected_atom != index:
                        atom_colors.append(atom_colors_dict[atom.atom_type])
                    elif index == self.selected_atom:
                        atom_colors.append([57, 255, 20, 255])
                    else:
                        atom_colors.append(default_atom_map[atom.atom_type])

                self.renderer = self.viewer.add_renderer(MyBallAndStickRenderer,
                                                         scale * coordinates,
                                                         np.array(type_array),
                                                         np.array(atom_colors),
                                                         np.array([bond.get_indices() for bond in self.bonds]),
                                                         bond_colors,

                                                         np.array(interactions),
                                                         np.array(interaction_colors)
                                                         )
            except Exception as e:
                self.disable(e)
