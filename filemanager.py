__author__ = 'martin'

import numpy as np
import time
import os

class FileManager:
    def __init__(self):
        self.last_write = time.time()
        self.start_time = time.time()

    def write_numpy_to_xyz(self, filename, atom_coordinates, atom_types):
        self.last_write = time.time()
        print("Writing to: "+filename)
        file = open(filename, "w")
        xyz_string = self.numpy_to_xyz_string(filename,atom_coordinates, atom_types)
        file.write(xyz_string)
        file.close()

    def numpy_to_xyz_string(self, filename, atom_coordinates, atom_types):
        signal_info = atom_types
        xyz_string = str(len(atom_coordinates))+"\n"
        xyz_string += "Energy Minimised Cartesian Coordinates\n"
        for i,atom in enumerate(atom_coordinates):
            xyz_string += (signal_info[i]+" "+str(atom_coordinates[i][0]/10)+" "+str(atom_coordinates[i][1]/10)+" "+str(atom_coordinates[i][2]/10)+"\n")
        return xyz_string

    def read_numpy_from_xyz(self, filename):
        xyzfile = open(filename, "r")
        xyz_data = xyzfile.readlines()
        xyzfile.close()
        xyz_data = xyz_data[2:]
        xyz_data = [line[2:] for line in xyz_data]
        xyz_data = [line.replace(" ",",") for line in xyz_data]
        xyz_data = ["["+line+"]" for line in xyz_data]
        xyz_data = [eval(line) for line in xyz_data]
        return np.array(xyz_data)

    def convert_xyz_to_mol(self, file_in):
        file_out = file_in[:-3]+"mol"
        os.system("obabel %s -O %s"%(file_in, file_out))

    def time_since_last_write(self):
        return time.time()-self.last_write

    def get_running_time(self):
        return time.time()-self.start_time


import main
if __name__ == "__main__":
    structure_solver = main.StructureMinimiser()
    interaction_manager = structure_solver.get_interaction_manager()
    file_manager = FileManager(interaction_manager)
    a = file_manager.read_numpy_from_xyz('Run2/solution1.xyz')