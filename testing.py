__author__ = 'martin'
from main_testing import get_interaction_manager
import numpy as np
from main_testing import magnitude
interaction_manager = get_interaction_manager()

def get_response_value(atom_coordinates, interaction_manager):
    """
    Calculate the response value
    """
    atom_coordinates = atom_coordinates.reshape((interaction_manager.number_signals,3))
    response = 0
    for i, v1 in enumerate(atom_coordinates):
        for j, v2 in enumerate(atom_coordinates):
            if j < i:
                response += interaction_manager.interaction_response(i, j, magnitude(v2-v1))
    return response

def get_coordinates_from_file():
    f = open('t2hexenal.np')
    s = f.read()
    f.close()
    s = eval(s)
    return np.array(s,dtype=float)

atom_coordinates = get_coordinates_from_file()
response_value = get_response_value(atom_coordinates,interaction_manager)
print(response_value)