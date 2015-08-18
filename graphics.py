from chemlab.graphics import QtViewer
from chemlab.graphics.renderers import BallAndStickRenderer
import numpy as np

class MolecularGraphics:
    def __init__(self, coordinates, type_array, bonds):
        self.viewer = QtViewer()
        self.active=True
        try:
            if len(coordinates)!=len(type_array):
                raise Exception

            self.type_array = np.array(type_array)
            self.coordinates = coordinates
            self.bonds = np.array(bonds)

            self.renderer = self.viewer.add_renderer(BallAndStickRenderer, 2*coordinates, self.type_array, self.bonds)
        except:
            self.active = False

    def update_coordinates(self, coordinates):
        if self.active:
            try:
                self.coordinates = coordinates
                self.renderer.update_positions(2*coordinates)
                self.viewer.widget.update()
            except:
                self.active=False


    def run(self, iterfunction):
        if self.active:
            try:
                self.viewer.schedule(iterfunction)
                self.viewer.run()
            except:
                self.active=False
        if not self.active:
            while True:
                iterfunction()

    def update(self, coordinates, type_array, bonds):
        if self.active:
            try:
                self.type_array = np.array(type_array)
                self.bonds = np.array(bonds)
                self.viewer.clear()
                self.renderer = self.viewer.add_renderer(BallAndStickRenderer, 2*self.coordinates, self.type_array, self.bonds)

            except:
                self.active=False




