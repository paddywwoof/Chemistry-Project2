__author__ = 'martin'

import numpy as np

class FileManager:
    def __init__(self, interaction_manager, resources_folder="resources/"):
        self.resources_folder = resources_folder
        self.interaction_manager = interaction_manager


    def write_numpy_to_xyz(self,filename, best_atom_coordinates):
        print("Writing to: "+self.resources_folder+filename)

        signal_info = self.interaction_manager.atom_types

        file = open(self.resources_folder+filename, "w")
        file.write(str(self.interaction_manager.number_signals)+"\n")
        file.write("Optimised cartesian coordinates\n")
        for i in range(self.interaction_manager.number_signals):
            file.write(signal_info[i]+" "+str(best_atom_coordinates[i][0]/10)+" "+str(best_atom_coordinates[i][1]/10)+" "+str(best_atom_coordinates[i][2]/10)+"\n")
        file.close()

    def read_xyz_to_numpy(self,filename):
        pass
        xyzfile = open(self.resources_folder+filename, "r")
        xyz_data = xyzfile.readlines()
        xyzfile.close()
        xyz_data = xyz_data[2:]
        xyz_data = [line[2:] for line in xyz_data]
        xyz_data = [line.replace(" ",",") for line in xyz_data]
        xyz_data = ["["+line+"]" for line in xyz_data]
        xyz_data = [eval(line) for line in xyz_data]
        return np.array(xyz_data)