from interactionmanager.interactionmanager import InteractionManager
from structureminimiser import StructureMinimiser
from signalmanager import OneDSignalManager, TwoDSignalManager
from filemanager import FileManager
import filemanager

def get_twod_signal_manager():
    oned_signal_manager = OneDSignalManager()
    oned_signal_manager.add_nmr_signals('resources/nmr/oned/hydrogen_integration_data.txt', "H")
    oned_signal_manager.add_nmr_signals('resources/nmr/oned/carbon_integration_data.txt', "C")
    twod_signal_manager = TwoDSignalManager(oned_signal_manager)
    twod_signal_manager.add_nmr_signals("COSY", 'resources/nmr/twod/cosy/cosy_peak_data.txt')
    twod_signal_manager.add_nmr_signals("HSQC", 'resources/nmr/twod/hsqc/hsqc_peak_data.txt')
    twod_signal_manager.add_nmr_signals("HMBC", 'resources/nmr/twod/hmbc/hmbc_peak_data.txt')
    return twod_signal_manager

def get_interaction_manager(interaction_matrix, type_array, shift_data):
    interaction_manager = InteractionManager(1001, interaction_matrix, type_array, shift_data)
    interaction_manager.add_default_interaction(0, "Default Repulsive", repulsive_amplitude=0.8, repulsive_time_constant=0.5, depth=10)
    interaction_manager.add_new_interaction(index=1, interaction_name="COSY 3-Bond H-H          ", repulsive_amplitude=0.8, repulsive_time_constant=6.28, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=2, interaction_name="HSQC 1-Bond H-C          ", repulsive_amplitude=0.8, repulsive_time_constant=2.26, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=3, interaction_name="HMBC 2/3 Bond H-C        ", repulsive_amplitude=0.8, repulsive_time_constant=6.50, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=4, interaction_name="INAD Single 1-Bond C-C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.49, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=5, interaction_name="INAD Double 1-Bond C=C   ", repulsive_amplitude=0.8, repulsive_time_constant=3.15, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.add_new_interaction(index=6, interaction_name="NOESY Normalised Function", repulsive_amplitude=0.8, repulsive_time_constant=2.03, depth=3, attractive_amplitude=0.6, attractive_time_constant=200, power=3)
    interaction_manager.plot_all_interactions()
    return interaction_manager

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
    test_main()
