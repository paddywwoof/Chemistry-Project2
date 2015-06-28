__author__ = 'martin'
import math
import random
import time

import scipy.optimize as opt

from interactions import InteractionManager
from get_interactions import get_hexenal_interactions

last_write = time.time()


def get_interaction_manager():
    interaction_manager = InteractionManager(number_carbon_signals=6, number_hydrogen_signals=10,
                                             get_interactions=get_hexenal_interactions, axis_width=1001)
    interaction_manager.add_default_interaction(0, repulsive_amplitude=0.8, repulsive_time_constant=5)
    interaction_manager.add_new_interaction(1, "COSY 3-Bond H-H   ", repulsive_amplitude=0.8, repulsive_time_constant=8.2, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(2, "HSQC 1-Bond H-C   ", repulsive_amplitude=0.8, repulsive_time_constant=2.3, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(3, "HMBC 2/3 Bond H-C ", repulsive_amplitude=0.8, repulsive_time_constant=6.5, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(4, "INAD 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.4, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    return interaction_manager


def get_response_value(atom_coordinates, interaction_manager):
    global last_write
    """
    Calculate the response value
    """
    atom_coordinates = atom_coordinates.reshape((interaction_manager.number_signals, 3))
    response = 0
    for i, v1 in enumerate(atom_coordinates):
        for j, v2 in enumerate(atom_coordinates):
            if j < i:
                response += interaction_manager.interaction_response(i, j, magnitude(v2-v1))

    if time.time()-last_write > 0.5:
        print(response)
        write_coordinates_to_xyz("resources/tempfile.xyz", interaction_manager, atom_coordinates)
        last_write = time.time()
    return response


def magnitude(vec):
    distance = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    if distance > 1000:
        distance = 1000
    return round(distance)


def main():
    start_time = time.time()
    interaction_manager = get_interaction_manager()
    interaction_manager.plot_all_interactions()
    best_atom_coordinates = interaction_manager.get_initial_coordinates()

    best_response_value = get_response_value(best_atom_coordinates, interaction_manager)

    number_attempts = 5  # Hardcoded

    for attempt_number in range(number_attempts):
        best_atom_coordinates = simulated_annealing(best_atom_coordinates, interaction_manager, get_response_value)

    best_atom_coordinates = best_atom_coordinates.reshape((interaction_manager.number_signals, 3))
    write_coordinates_to_xyz("resources/tempfile.xyz", interaction_manager, best_atom_coordinates)
    write_coordinates_to_xyz("resources/solution.xyz", interaction_manager, best_atom_coordinates)
    running_time = time.time()-start_time
    print("Running time: ", running_time, "s")
    return best_response_value, best_atom_coordinates


def write_coordinates_to_xyz(filename, interaction_manager, best_atom_coordinates):
    print("Writing to: "+filename)
    signal_info = ["H"]*interaction_manager.number_hydrogen_signals + ["C"]*interaction_manager.number_carbon_signals
    file = open(filename, "w")
    file.write(str(interaction_manager.number_signals)+"\n")
    file.write("Optimised cartesian coordinates\n")
    for i in range(interaction_manager.number_signals):
        file.write(signal_info[i]+" "+str(best_atom_coordinates[i][0]/10)+" "+str(best_atom_coordinates[i][1]/10)+" "+str(best_atom_coordinates[i][2]/10)+"\n")
    file.close()


if __name__ == "__main__":
    response_value, coordinates = main()




