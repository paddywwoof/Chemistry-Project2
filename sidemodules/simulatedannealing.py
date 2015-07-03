import random
from math import exp
step_size = 0.5

def simulated_annealing(atom_coordinates, interaction_manager, get_response_value):

    t = 30
    b = 1.001

    while t>0.00001:

        t /= b

        for atom in atom_coordinates:
            pass

        response_value = get_response_value(atom_coordinates,interaction_manager)
        index = random.randrange(0, len(atom_coordinates))
        atom = atom_coordinates[index]
        delta = [0,0,0]
        delta[random.randrange(0,3)] = random.uniform(0, step_size)
        atom += delta
        new_response_value = get_response_value(atom_coordinates, interaction_manager)

        deltaE = new_response_value - response_value

        if deltaE < 0:
            continue
        elif deltaE >=0 and t>0:
            if random.uniform(0, 1) < exp(-deltaE*100.0/t):
                continue
        atom -= delta
    return atom_coordinates












