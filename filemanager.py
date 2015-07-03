__author__ = 'martin'

import numpy as np
import time

class FileManager:
    def __init__(self, interaction_manager, resources_folder="resources/"):
        self.resources_folder = resources_folder
        self.interaction_manager = interaction_manager
        self.last_write = time.time()
        self.start_time = time.time()

    def write_numpy_to_xyz(self, filename, best_atom_coordinates):
        self.last_write = time.time()
        print("Writing to: "+self.resources_folder+filename)
        file = open(self.resources_folder+filename, "w")
        xyz_string = self.numpy_to_xyz_string(filename,best_atom_coordinates)
        file.write(xyz_string)
        file.close()

    def numpy_to_xyz_string(self, filename, best_atom_coordinates):
        signal_info = self.interaction_manager.atom_types
        xyz_string = ""
        xyz_string += (str(self.interaction_manager.number_signals)+"\n")
        xyz_string += ("Optimised cartesian coordinates\n")
        for i in range(self.interaction_manager.number_signals):
            xyz_string += (signal_info[i]+" "+str(best_atom_coordinates[i][0]/10)+" "+str(best_atom_coordinates[i][1]/10)+" "+str(best_atom_coordinates[i][2]/10)+"\n")
        return xyz_string

    def read_xyz_to_numpy(self,filename):
        xyzfile = open(self.resources_folder+filename, "r")
        xyz_data = xyzfile.readlines()
        xyzfile.close()
        xyz_data = xyz_data[2:]
        xyz_data = [line[2:] for line in xyz_data]
        xyz_data = [line.replace(" ",",") for line in xyz_data]
        xyz_data = ["["+line+"]" for line in xyz_data]
        xyz_data = [eval(line) for line in xyz_data]
        return np.array(xyz_data)

    def time_since_last_write(self):
        return time.time()-self.last_write

    def get_running_time(self):
        return time.time()-self.start_time


import main
structure_solver = main.StructureMinimiser()
interaction_manager = structure_solver.get_interaction_manager()
file_manager = FileManager(interaction_manager)
a = file_manager.read_xyz_to_numpy('Run2/solution1.xyz')
print(a)