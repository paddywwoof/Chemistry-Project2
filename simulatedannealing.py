import random

step_size = 1

def simulated_annealing(atom_coordinates, interaction_manager, get_response_value):

    t = 100
    b = 1.01

    for x in range(0, 10000):
        t /= b

        response_value = get_response_value(atom_coordinates,interaction_manager)
        index = random.randrange(0, len(atom_coordinates))
        atom = atom_coordinates[index]
        delta = [random.uniform(0, step_size) for i in range(0, 3)]
        atom += delta
        new_response_value = get_response_value(atom_coordinates, interaction_manager)

        deltaE = new_response_value - response_value

        if deltaE<0:
            continue
        elif T>0:
            if random.uniform(0,1)<exp(deltaE*100.0/t):
                continue
            else:
                atom -= delta
    return atom_coordinates












