from main_run import get_interaction_manager,get_twod_signal_manager
from filemanager import FileManager
from chemlab.graphics import QtViewer
import numpy as np
import random
from chemlab.graphics.transformations import rotation_matrix
from chemlab.graphics.renderers import BallAndStickRenderer
from math import exp
import os

minima = {0: 1.0, 1: 0.243, 2: 0.109, 3: 0.25, 4: 0.154, 5: 0.142, 6: 0.3}

class ChemLabMinimiser:
    def __init__(self):
        self.interaction_manager = get_interaction_manager(*get_twod_signal_manager().get_interaction_data())

        #TODO : REMOVE WEIRD POINTS
        for x in [[19,6],[19,7],[21,6]]:
            print(self.interaction_manager.interaction_matrix[x])
            self.interaction_manager.interaction_matrix[x[0]][x[1]] = 0
            self.interaction_manager.interaction_matrix[x[1]][x[0]] = 0

        self.interaction_matrix = self.interaction_manager.interaction_matrix
        self.interaction_matrix[self.interaction_matrix==9]=0
        self.type_array = self.interaction_manager.type_array
        self.shift_data = self.interaction_manager.shift_data
        self.number_atoms = self.interaction_manager.number_atoms

        self.viewer = QtViewer()

        ac = np.random.rand(self.number_atoms, 3)*3
        self.coordinate_manager = CoordinateManager(ac, self.interaction_manager)

        self.renderer = self.viewer.add_renderer(BallAndStickRenderer, self.coordinate_manager.get_coordinates(), self.type_array, np.array([]))


    def main(self):
        self.viewer.schedule(self.iteration)
        self.viewer.run()

    def iteration(self):
        self.force_step()

    def force_step(self):
        atom_coordinates = self.coordinate_manager.get_coordinates()
        for index, atom in enumerate(atom_coordinates):
            atom_coordinates[index] += self.get_force_vector(index, atom, atom_coordinates)
            self.renderer.update_positions(atom_coordinates)
        self.coordinate_manager.update(atom_coordinates)

    def get_force_vector(self, index, atom, atom_coordinates):
        vec = np.zeros(3)
        for i2, atom2 in enumerate(atom_coordinates):


            distance = np.linalg.norm(atom2-atom)
            if distance == 0 or atom is atom2:
                continue
            #input((atom,atom2))



            subvec = 0.01*(atom2-atom)/distance

            minimum = minima[self.interaction_matrix[index][i2]]

            try:
                delta = abs(minimum-distance)
            except Exception as e:
                pass




            try:
                #print(delta, minimum, distance)
                subvec *= (delta*delta)
                vec+=subvec
            except:
                pass




        return vec

class CoordinateManager:
    def __init__(self, atom_coordinates, interaction_manager):
        self.atom_coordinates = atom_coordinates

    def update(self, atom_coordinates):
        for x in [0, 1, 2]:
            atom_coordinates[:, x] -= min(atom_coordinates[:, x])
        self.atom_coordinates = atom_coordinates

    def get_coordinates(self):
        return np.copy(self.atom_coordinates)


def main():
    clm = ChemLabMinimiser()
    clm.main()
    input("Press Any Key to Quit")
    raise SystemExit

if __name__ == "__main__":
    main()
