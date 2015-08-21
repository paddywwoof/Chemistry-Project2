from chemlab.graphics.renderers.base import AbstractRenderer
from chemlab.graphics.renderers.atom import AtomRenderer
from chemlab.graphics.renderers.line import LineRenderer

from .mybondrenderer import MyBondRenderer


from chemlab.db import ChemlabDB
cdb = ChemlabDB()

class MyBallAndStickRenderer(AbstractRenderer):
    '''Render a ball and stick representation of a series of
    coordinates and bonds.
    
    .. image:: /_static/ballandstick_renderer.png
    
    **Parameters**
    
    widget:
        The parent QChemlabWidget
    r_array: np.ndarray((NATOMS, 3), dtype=float)
        The coordinate array
    type_array: np.ndarray((NATOMS, 3), dtype=object)
        An array containing all the atomic symbols like `Ar`, `H`, `O`.
        If the atomic type is unknown, use the `Xx` symbol.
    bonds: np.ndarray((NBONDS, 2), dtype=int)
        An array of integer pairs that represent the bonds.


    '''
    def __init__(self, widget, r_array, type_array, bonds, bond_colors, interactions, interaction_colors, shading='phong'):
        super(MyBallAndStickRenderer, self).__init__(widget)
        vdw_dict = cdb.get("data", 'vdwdict')
        
        scale = 0.3
        for k in vdw_dict:
            vdw_dict[k] = vdw_dict[k] * scale
        
        self.has_bonds = len(bonds) > 0
        
        self.ar = AtomRenderer(widget, r_array, type_array, radii_map=vdw_dict, shading=shading)

        if self.has_bonds:
            self.br = MyBondRenderer(widget, bonds, bond_colors, r_array, style='impostors', shading=shading)
        self.has_interactions = len(interactions) > 0

        if self.has_interactions:
            if len(interactions) != len(interaction_colors):
                raise Exception("Should be 1-1 correspondence between interaction and color")

            self.interaction_renderer = MyBondRenderer(widget, interactions, interaction_colors, r_array, style='impostors', shading=shading, radius=0.005)

    def draw(self):
        self.ar.draw()

        if self.has_bonds:
            self.br.draw()

        if self.has_interactions:
            self.interaction_renderer.draw()

    def update_positions(self, r_array):
        """Update the coordinate array r_array"""
        self.ar.update_positions(r_array)
        
        if self.has_bonds:
            self.br.update_positions(r_array)

        if self.has_interactions:
            self.interaction_renderer.update_positions(r_array)
