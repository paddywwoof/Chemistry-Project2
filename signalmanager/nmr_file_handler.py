from . import OneDSignalManager, TwoDSignalManager
import zipfile

def get_twod_signal_manager(path):
    nmr_file = zipfile.ZipFile(path, "r")

    nmr_string_dict = {}
    for nmr_type in ["PROTON", "CARBON", "COSY", "HMBC", "HSQC", "NOESY"]:
        file_name_list = [x for x in nmr_file.namelist() if x.split(".")[0].lower() == nmr_type.lower()]
        if len(file_name_list) > 0:
            file_name = file_name_list[0]
            nmr_string = nmr_file.read(file_name).decode(encoding='UTF-8')
            nmr_string_dict[nmr_type] = nmr_string
        else:
            nmr_string_dict[nmr_type] = ""
    nmr_file.close()

    oned_signal_manager = OneDSignalManager()
    for oned in [("PROTON", "H"), ("CARBON", "C")]:
        oned_signal_manager.add_nmr_signals(nmr_string_dict[oned[0]], oned[1])

    twod_signal_manager = TwoDSignalManager(oned_signal_manager)
    for twod in ["COSY", "HSQC", "HMBC", "NOESY"]:
        twod_signal_manager.add_nmr_signals(twod, nmr_string_dict[twod])

    return twod_signal_manager
