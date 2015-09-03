__author__ = 'martin'

import time
import os
import random
import numpy as np

class FileManager:
    def __init__(self):
        self.last_write = time.time()
        self.start_time = time.time()

    def write_numpy_to_xyz(self, filename, atom_coordinates, type_array):
        print("Writing to %s" % (os.path.abspath(filename)))
        self.last_write = time.time()
        xyz_string = self.convert_numpy_to_xyz_string(atom_coordinates, type_array)
        writefile(filename, xyz_string)

    def convert_numpy_to_xyz_string(self, atom_coordinates, type_array):
        signal_info = type_array
        xyz_string = str(len(atom_coordinates))+"\n"
        xyz_string += "Energy Minimised Cartesian Coordinates\n"
        for i, atom in enumerate(atom_coordinates):
            xyz_string += (signal_info[i]+" "+str(atom_coordinates[i][0]/10)+" "+str(atom_coordinates[i][1]/10)+" "+str(atom_coordinates[i][2]/10)+"\n")
        return xyz_string

    def read_numpy_from_xyz(self, filename):
        xyz_string = readfile(filename)
        xyz_data = xyz_string.splitlines()
        xyz_data = xyz_data[2:]
        for i, line in enumerate(xyz_data):
            while xyz_data[i][0] == " ":
                xyz_data[i] = xyz_data[i][1:]

            while "  " in xyz_data[i]:
                xyz_data[i] = xyz_data[i].replace("  ", " ")
            xyz_data[i] = xyz_data[i][2:]

        xyz_data = [line.replace(" ", ",") for line in xyz_data]

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
        type_array = []
        shift_values = []
        for atom in interaction_list:
            type_array.append(atom[0])
            interaction_matrix.append(eval(atom[1].replace(" ", ",")))
            shift_values.append(float(atom[2]))
        interaction_matrix = np.array(interaction_matrix)
        return interaction_matrix, type_array, shift_values

    def write_numpy_to_mol(self,filename, bonds, atoms, atom_coordinates = None):
        mol_file_string = self.convert_numpy_to_mol_string(bonds, atoms, atom_coordinates)
        file = open(filename, "w")
        file.write(mol_file_string)
        file.close()

    def convert_numpy_to_mol_string(self, bonds, atoms, atom_coordinates=None):
        header = """Molecule Name \n     Additional Information\n\n %s %s  0  0  0  0  0  0  0  0999 V2000\n""" % (len(atoms), len(bonds))
        footer = "M  END"
        line = "    %s    %s    %s %s "+"  0"*12
        mol_file_string = ""
        mol_file_string += header
        if atom_coordinates is None:
            for index, atom in enumerate(atoms):
                i1 = str(random.uniform(0, 1))[:6]
                i2 = str(random.uniform(0, 1))[:6]
                i3 = str(random.uniform(0, 1))[:6]
                mol_file_string += line % (i1, i2, i3, atom.atom_type)+"\n"
        else:
            for index, atom in enumerate(atoms):
                i1 = str(10*atom_coordinates[index][0])[:6] + "0" * max(0, 6-len(str(atom_coordinates[index][0])[:6]))
                i2 = str(10*atom_coordinates[index][1])[:6] + "0" * max(0, 6-len(str(atom_coordinates[index][1])[:6]))
                i3 = str(10*atom_coordinates[index][2])[:6] + "0" * max(0, 6-len(str(atom_coordinates[index][2])[:6]))
                mol_file_string += line % (i1, i2, i3, atom.atom_type)+"\n"
        bond_line = " %s %s  %s  0  0  0  0"
        for index, bond in enumerate(bonds):
            pass
            a1 = " "*(2-len(str(bond[0]+1)))+str(bond[0]+1)
            a2 = " "*(2-len(str(bond[1]+1)))+str(bond[1]+1)
            a3 = bond.bond_order

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