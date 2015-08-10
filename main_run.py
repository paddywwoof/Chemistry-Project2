from interactionmanager.interactionmanager import InteractionManager
from structureminimiser import StructureMinimiser
from signalmanager import OneDSignalManager, TwoDSignalManager
from filemanager import FileManager
import filemanager


def main():
    twod_signal_manager = get_twod_signal_manager()
    interaction_manager = get_interaction_manager(*twod_signal_manager.get_interaction_data())


    for i, vec in enumerate(interaction_manager.interaction_matrix):
        for j, vec2 in enumerate(vec):
            if interaction_manager.interaction_matrix[i][j]==4:
                print([i+1,j+1])

    interaction_manager.print_matrix()
    structure_minimiser = StructureMinimiser(interaction_manager)
    file_manager = FileManager()

    mol_string = file_manager.convert_numpy_to_mol_string(interaction_manager)
    filemanager.writefile('outfile.mol', mol_string)

    try:
        structure_minimiser.minimise_response()
    except KeyboardInterrupt:
        pass
    structure_minimiser.write_solution("solution.xyz")
    interaction_manager.print_matrix()


def test_main():
    interaction_manager = InteractionManager(1001, *get_twod_signal_manager().get_interaction_data())
    for x in range(1, 2000):
        print("RTC = %s"%(x/250))
        interaction_manager.add_new_interaction(index=1, interaction_name="COSY 3-Bond H-H   ", repulsive_amplitude=0.8, repulsive_time_constant=x/250, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    for x in range(1, 2000):
        print("RTC = %s"%(x/250))
        interaction_manager.add_new_interaction(index=2, interaction_name="HSQC 1-Bond H-C   ", repulsive_amplitude=0.8, repulsive_time_constant=x/250, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    for x in range(1, 2000):
        print("RTC = %s"%(x/250))
        interaction_manager.add_new_interaction(index=3, interaction_name="HMBC 2/3 Bond H-C ", repulsive_amplitude=0.8, repulsive_time_constant=x/250, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    for x in range(1, 2000):
        print("RTC = %s"%(x/250))
        interaction_manager.add_new_interaction(index=4, interaction_name="INAD 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=x/250, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)

    return interaction_manager


if __name__ == "__main__":
    main()
