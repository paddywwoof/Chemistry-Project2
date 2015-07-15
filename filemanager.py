__author__ = 'martin'

import time
import os
import random
import numpy as np

class FileManager:
    def __init__(self):
        self.last_write = time.time()
        self.start_time = time.time()

    def write_numpy_to_xyz(self, filename, atom_coordinates, atom_types):
        self.last_write = time.time()
        xyz_string = self.convert_numpy_to_xyz_string(atom_coordinates, atom_types)
        writefile(filename, xyz_string)

    def convert_numpy_to_xyz_string(self, atom_coordinates, atom_types):
        signal_info = atom_types
        xyz_string = str(len(atom_coordinates))+"\n"
        xyz_string += "Energy Minimised Cartesian Coordinates\n"
        for i,atom in enumerate(atom_coordinates):
            xyz_string += (signal_info[i]+" "+str(atom_coordinates[i][0]/10)+" "+str(atom_coordinates[i][1]/10)+" "+str(atom_coordinates[i][2]/10)+"\n")
        return xyz_string

    def read_numpy_from_xyz(self, filename):
        xyz_string = readfile(filename)
        xyz_data = xyz_string.splitlines()
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

    def reset(self):
        self.start_time = time.time()
        self.last_write = time.time()

    def load_interaction_data(self, filename):
        interaction_string = readfile(filename)
        interaction_list = [line.split(":") for line in interaction_string.splitlines()]
        interaction_matrix = []
        atom_types = []
        shift_values = []
        for atom in interaction_list:
            atom_types.append(atom[0])
            interaction_matrix.append(eval(atom[1].replace(" ", ",")))
            shift_values.append(float(atom[2]))
        interaction_matrix = np.array(interaction_matrix)
        return interaction_matrix, atom_types, shift_values

    def write_numpy_to_mol(self,filename, interaction_manager):
        mol_file_string = self.convert_numpy_to_mol_string(interaction_manager)
        file = open(filename, "w")
        file.write(mol_file_string)
        file.close()

    def convert_numpy_to_mol_string(self, interaction_manager):
        atom_types = interaction_manager.atom_types
        bonds = interaction_manager.bonds

        header = """Molecule Name \n     Additional Information\n\n     %s %s  0  0  0  0  0  0  0  0999 V2000\n    """ % (len(atom_types), len(bonds))
        footer = "M  END"
        line = "    %s    %s    %s %s "+"  0"*12
        mol_file_string = ""
        mol_file_string += header
        for atom in atom_types:
            i1 = str(random.uniform(0, 1))[:6]
            i2 = str(random.uniform(0, 1))[:6]
            i3 = str(random.uniform(0, 1))[:6]
            mol_file_string += line % (i1, i2, i3, atom)+"\n"
        bond_line = " %s %s  %s  0  0  0  0"
        for bond in bonds:
            pass
            a1 = " "*(2-len(str(bond[0])))+str(bond[0])
            a2 = " "*(2-len(str(bond[1])))+str(bond[1])
            a3 = bond[2]
            b = bond_line % (a1, a2, a3)
            mol_file_string += (b+"\n")
        mol_file_string += footer
        return mol_file_string

def readfile(path):
    current_path = os.path.dirname(__file__) + "/"
    openfile = open(current_path+path, "r")
    file_string = openfile.read()
    openfile.close()
    return file_string

def writefile(path, string):
    file = open(path, "w")
    file.write(string)
    file.close()