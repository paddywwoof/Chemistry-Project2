from chemlab.mviewer.qtmolecularviewer import QtMolecularViewer
from chemlab.mviewer.representations.ballandstick import BallAndStickRepresentation
from chemlab.core import System, Atom, Molecule
import numpy as np

class MyBallAndStickRepresentation(BallAndStickRepresentation):
    def on_atom_selection_changed(self):
        sel = self.selection_state['atoms'].indices
        cols = np.array(self.atom_colors.array)
        cols[sel, -1] = 50
        self.atom_renderer.update_colors(cols)
        self.viewer.update()


a1 = Atom('Ar', [0.0, 0.0, 0.0])
a2 = Atom('N', [0.0, 0.5, 1.0])
mol = Molecule([a1,a2])
system = System([mol])

viewer = QtMolecularViewer()
viewer.add_representation(MyBallAndStickRepresentation, system)

def print_indices():
    indices = viewer.representation.selection_state['atoms'].indices
    if len(indices)>0:
        print(indices)

viewer.schedule(print_indices)
viewer.run()

