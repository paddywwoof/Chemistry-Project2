__author__ = 'martin'
from math import sqrt
from scipy.optimize import basinhopping
from interactionmanager import InteractionManager
from filemanager import FileManager


class StructureMinimiser:
    def __init__(self):
        self.interaction_manager = self.get_interaction_manager()
        self.interaction_manager.plot_all_interactions()

        self.file_manager = FileManager(self.interaction_manager)

        self.atom_coordinates = self.interaction_manager.get_initial_coordinates()
        self.response_value = 0

    def main(self):
        import os
        for run in range(3,10):
            print("Starting Run: ", run)
            new_dir = "Run"+str(run)
            try:
                os.mkdir("resources/"+new_dir)
            except FileExistsError:
                pass
            self.minimise_response()
            self.write_solution(new_dir+"/")
            print("Running time: ", self.file_manager.get_running_time(), "s")
        return self.response_value, self.atom_coordinates

    def write_solution(self, directory=""):
        for x in [1.15, 1.1, 1.05, 1, 0.95, 0.85, 0.80]:
            self.file_manager.write_numpy_to_xyz(directory+"solution%s.xyz" % x, self.atom_coordinates * x)

    def minimise_response(self):
        minimisation_solution = basinhopping(self.calculate_response_value, x0=self.atom_coordinates, niter=50, minimizer_kwargs={"method": "Nelder-Mead"}, T=30, stepsize=0.5)
        self.response_value = minimisation_solution.fun
        self.atom_coordinates = minimisation_solution.x.reshape((self.interaction_manager.number_signals, 3))

    def get_interaction_manager(self):
        interaction_manager = InteractionManager(axis_width=1001, interaction_filename="interaction_cache/cucumber.np")
        interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=5)
        interaction_manager.add_new_interaction(index=1, interaction_name="COSY 3-Bond H-H   ", repulsive_amplitude=0.8, repulsive_time_constant=8.2, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
        interaction_manager.add_new_interaction(index=2, interaction_name="HSQC 1-Bond H-C   ", repulsive_amplitude=0.8, repulsive_time_constant=2.3, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
        interaction_manager.add_new_interaction(index=3, interaction_name="HMBC 2/3 Bond H-C ", repulsive_amplitude=0.8, repulsive_time_constant=6.5, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
        interaction_manager.add_new_interaction(index=4, interaction_name="INAD 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.4, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
        return interaction_manager

    def calculate_response_value(self, atom_coordinates, write_out=True):

        atom_coordinates = atom_coordinates.reshape((self.interaction_manager.number_signals, 3))
        response = 0
        for i, v1 in enumerate(atom_coordinates):
            for j, v2 in enumerate(atom_coordinates):
                if j < i:
                    response += self.interaction_manager.interaction_response(i, j, self.calculate_magnitude(v2 - v1))

        if self.file_manager.time_since_last_write() > 0.5 and write_out:
            print("Response Value: ", response)
            self.file_manager.write_numpy_to_xyz("tempfile.xyz", atom_coordinates)
        return response

    def calculate_magnitude(self, vec):
        if vec.tolist() == [0, 0, 0]:
            return 0
        try:
            distance = sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)
        except Exception as e:
            print("Vector:" + str(vec))
            raise e
        return round(distance)


if __name__ == "__main__":
    try:
        response_value, coordinates = StructureMinimiser().main()
    except KeyboardInterrupt:
        print("\n \n \n Program Terminated By User")
