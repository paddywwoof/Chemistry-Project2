import numpy as np


def get_hexenal_interactions():
    return get_interaction('hexenal.np')

def get_cucumber_interactions():
    return get_interaction('cucumber.np')


def get_interaction(name):
    array_file = open("resources/"+name)
    array_string = array_file.read()
    array_string = array_string.replace(" ", ", ")
    array_file.close()
    interaction_matrix = eval(array_string)
    interaction_matrix = np.array(interaction_matrix)
    interaction_matrix = interaction_matrix.transpose()
    return interaction_matrix
