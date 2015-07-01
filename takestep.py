import numpy as np
import random
class SearchCube:
    def __init__(self, width, length, height, calculate_response, interaction_manager):
        self.width=width
        self.length = length
        self.height = height
        self.calculate_response = calculate_response
        self.interaction_manager = interaction_manager

    def __call__(self, coords):
        """take a random step but ensure the new position is within the bounds"""
        coords = coords.reshape((self.interaction_manager.number_signals, 3))

        min_coords = coords

        min_response = self.calculate_response(coords, self.calculate_response)

        atom_index = random.randrange(0,len(coords))

        for x in range(-width,width):
            for y in range(-length,length):
                for z in range(-height,height):
                    coords[atom_index]=[x,y,z]
                    response = self.calculate_response(coords, self.calculate_response)
                    if response < min_response:
                        min_response = response
                        min_coords = np.copy(coords)
        return min_coords