__author__ = 'martin'
import math
import random
import time
import scipy.optimize as opt

from interactions import InteractionManager
from get_interactions import get_cucumber_interactions as interaction_getter
import numpy as np
"""
Interaction_manager parameters
"""

last_write = time.time()


def get_interaction_manager():
    interaction_manager = InteractionManager(get_interactions=interaction_getter, axis_width=1001)

    interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=5, depth=1.0)
    interaction_manager.add_new_interaction(1, "COSY 3-Bond H-H   ", repulsive_amplitude=0.8, repulsive_time_constant=8.2, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(2, "HSQC 1-Bond H-C   ", repulsive_amplitude=0.8, repulsive_time_constant=2.3, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(3, "HMBC 2/3 Bond H-C ", repulsive_amplitude=0.8, repulsive_time_constant=6.5, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(4, "INAD 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.4, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)

    """
    interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=5, depth=1.0)
    interaction_manager.add_new_interaction(1, "COSY 3-Bond H-H   ", repulsive_amplitude=0.8, repulsive_time_constant=7.2, depth=1.0, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(2, "HSQC 1-Bond H-C   ", repulsive_amplitude=0.8, repulsive_time_constant=1.8, depth=1.0, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(3, "HMBC 2/3 Bond H-C ", repulsive_amplitude=0.8, repulsive_time_constant=5.5, depth=1.0, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(4, "INAD 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=2.7, depth=1.0, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    """
    """
    for x in range(1,500):
        y = round(x*0.1,3)
        interaction_manager.add_new_interaction(1, "COSY 3-Bond H-H "+str(y), repulsive_amplitude=0.8, repulsive_time_constant=y, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    for x in range(1,500):
        y = round(x*0.1,3)
        interaction_manager.add_new_interaction(2, "HSQC 1-Bond H-C "+str(y), repulsive_amplitude=0.8, repulsive_time_constant=y, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    for x in range(1,500):
        y = round(x*0.1,3)
        interaction_manager.add_new_interaction(3, "HMBC 2/3 Bond H-C "+str(y), repulsive_amplitude=0.8, repulsive_time_constant=y, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    for x in range(1,500):
        y = round(x*0.1,3)
        interaction_manager.add_new_interaction(4, "INAD 1-Bond C-C "+str(y), repulsive_amplitude=0.8, repulsive_time_constant=y, depth=1, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    """
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

    if time.time()-last_write > 0.1:
        print("Response Value: ",response)


        write_coordinates_to_xyz("resources/tempfile.xyz", interaction_manager, atom_coordinates)
        last_write = time.time()
    return response*3


def magnitude(vec):
    if vec.tolist() == [0,0,0]:
        return 0
    try:
        distance = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    except Exception as e:
        print("Vector:"+str(vec))
        raise e
    return round(distance)


def main():
    start_time = time.time()
    interaction_manager = get_interaction_manager()
    interaction_manager.plot_all_interactions()
    best_atom_coordinates = interaction_manager.get_initial_coordinates()

    for x in range(0, 50):
        if x%10==0:
            print("Iteration: ",x)
            print(best_atom_coordinates)
            write_coordinates_to_xyz("resources/tempfile.xyz", interaction_manager, best_atom_coordinates)
        for atom1_index, atom1 in enumerate(best_atom_coordinates):
            for atom2_index,atom2 in enumerate(best_atom_coordinates):

                if atom1.tolist() != atom2.tolist():
                    i_type = interaction_manager.interaction_matrix[atom1_index][atom2_index]
                    if i_type == 1:
                        if magnitude(atom1 - atom2) > 32*2:
                            best_atom_coordinates[atom2_index] = atom1 + (atom2-atom1)*0.33
                            best_atom_coordinates[atom1_index] = atom2 + (atom1-atom2)*0.33
                    if i_type == 2:
                        if magnitude(atom1 - atom2) > 12*2:
                            best_atom_coordinates[atom2_index] = atom1 + (atom2-atom1)*0.33
                            best_atom_coordinates[atom1_index] = atom2 + (atom1-atom2)*0.33
                    if i_type == 3:
                        if magnitude(atom1 - atom2) > 26*2:
                            best_atom_coordinates[atom2_index] = atom1 + (atom2-atom1)*0.33
                            best_atom_coordinates[atom1_index] = atom2 + (atom1-atom2)*0.33
                    if i_type == 4:
                        if magnitude(atom1 - atom2) > 16*2:
                            best_atom_coordinates[atom2_index] = atom1 + (atom2-atom1)*0.33
                            best_atom_coordinates[atom1_index] = atom2 + (atom1-atom2)*0.33

    write_coordinates_to_xyz("resources/tempfile.xyz", interaction_manager, best_atom_coordinates)
    input("")
    best_response_value = get_response_value(best_atom_coordinates, interaction_manager)

    number_attempts = 5  # Hardcoded

    for attempt_number in range(number_attempts):
        try:
            minimizer_kwargz = {
                "method": "Nelder-Mead",
                "args": (interaction_manager,)
                }

            """
            best_atom_coordinates = best_atom_coordinates.reshape((interaction_manager.number_signals, 3))

            for atom1_index,atom1 in enumerate(best_atom_coordinates):

                if max([magnitude(atom1-atom2) for atom2 in best_atom_coordinates if atom1 is not atom2])>50:
                    if atom1_index>0:
                        best_atom_coordinates[atom1_index] = best_atom_coordinates[atom1_index-1]
                    else:
                        best_atom_coordinates[atom1_index] = best_atom_coordinates[atom1_index+1]

            """
            """
            best_atom_coordinates*=1.01
            while get_response_value(best_atom_coordinates*0.99,interaction_manager) < get_response_value(best_atom_coordinates,interaction_manager):
                best_atom_coordinates = best_atom_coordinates * 0.99
                print(best_atom_coordinates)
            """


            out_min = opt.basinhopping(get_response_value, x0=best_atom_coordinates, niter=100, minimizer_kwargs=minimizer_kwargz, T=30, stepsize=0.5)

            attempt_atom_coordinates = out_min.x
            attempt_response_value = out_min.fun

            if attempt_response_value < best_response_value:
                best_atom_coordinates = attempt_atom_coordinates
                best_response_value = attempt_response_value
        except KeyboardInterrupt:
            break



    best_atom_coordinates = best_atom_coordinates.reshape((interaction_manager.number_signals, 3))
    write_coordinates_to_xyz("resources/tempfile.xyz", interaction_manager, best_atom_coordinates)
    write_coordinates_to_xyz("resources/solution.xyz", interaction_manager, best_atom_coordinates)
    running_time = time.time()-start_time
    print("Running time: ", running_time, "s")
    return best_response_value, best_atom_coordinates


def write_coordinates_to_xyz(filename, interaction_manager, best_atom_coordinates):
    print("Writing to: "+filename)

    signal_info = interaction_manager.atom_types

    file = open(filename, "w")
    file.write(str(interaction_manager.number_signals)+"\n")
    file.write("Optimised cartesian coordinates\n")
    for i in range(interaction_manager.number_signals):
        file.write(signal_info[i]+" "+str(best_atom_coordinates[i][0]/10)+" "+str(best_atom_coordinates[i][1]/10)+" "+str(best_atom_coordinates[i][2]/10)+"\n")
    file.close()


if __name__ == "__main__":
    response_value, coordinates = main()




