__author__ = 'martin'
import os
import datetime
from scipy.optimize import basinhopping
import numpy as np
from filemanager import FileManager


class StructureMinimiser:
    def __init__(self, interaction_manager):
        self.interaction_manager = interaction_manager
        self.interaction_manager.plot_all_interactions()
        self.file_manager = FileManager()
        self.best_atom_coordinates = self.interaction_manager.get_initial_coordinates()
        self.best_response_value = 0
        self.iterations = 0
        self.date = self.get_now()

    def reset(self):
        self.best_atom_coordinates = self.interaction_manager.get_initial_coordinates()
        self.best_response_value = 0
        self.iterations = 0
        self.file_manager.reset()

    def write_solution(self, filename, scale):
        self.file_manager.write_numpy_to_xyz("resources/" + filename, self.best_atom_coordinates * scale, self.interaction_manager.type_array)
        self.file_manager.convert_xyz_to_mol("resources/" + filename)

    def minimise_response(self):
        self.reset()
        basinhopping(self.calculate_response_value, x0=self.best_atom_coordinates, niter=1000, minimizer_kwargs={"method": "Nelder-Mead"}, T=30, stepsize=5)
        print("Running time: ", self.file_manager.get_running_time(), "s")

    def calculate_response_value(self, atom_coordinates, write_out=True):
        self.iterations += 1
        atom_coordinates = atom_coordinates.reshape((self.interaction_manager.number_signals, 3))
        response = 0
        for i, v1 in enumerate(atom_coordinates):
            for j, v2 in enumerate(atom_coordinates):
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, v2, v1)

        if response < self.best_response_value:
            self.best_response_value = response
            self.best_atom_coordinates = np.copy(atom_coordinates)

        if self.file_manager.time_since_last_write() > 0.5 and write_out:
            speed = self.iterations*1.0/self.file_manager.get_running_time()
            print("Response Value: ", response, " Iterations: ", self.iterations, "Speed: ", speed, end="\r")
            self.file_manager.write_numpy_to_xyz("output/tempfile.xyz", atom_coordinates, self.interaction_manager.type_array)
        return response

    def get_now(self):
        date = str(datetime.datetime.now())[:-7]
        date = date.replace(" ", "~")
        date = date.replace(":", "-")
        return date


def minimisation_test(interaction_manager):
    structure_minimiser = StructureMinimiser(interaction_manager)
    try:
        os.mkdir("output")
    except FileExistsError:
        pass
    try:
        os.mkdir("output/"+structure_minimiser.date)
    except FileExistsError:
        pass
    for run in range(1, 11):
        print("Starting Run: ", run)
        new_dir = "output/"+structure_minimiser.date+"/Run"+str(run)+"/"
        try:
            os.mkdir(new_dir)
        except FileExistsError:
            pass
        structure_minimiser.minimise_response()
        for scale in [1.15, 1.1, 1.05, 1, 0.95, 0.85, 0.80]:
            structure_minimiser.write_solution("solution%s.xyz" % scale)

if __name__ == "__main__":
    minimisation_test()
