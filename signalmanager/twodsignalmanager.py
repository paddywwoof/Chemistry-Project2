import numpy as np
import matplotlib.pyplot as plt
from itertools import product as cartesian
from .signalparsers import parse_table

np.set_printoptions(suppress=True)

class InteractionValues:
    DEFAULT = 0
    COSY = 1  # Through Space
    HSQC = 2  # Bond Type
    HMBC = 3  # Through Space
    INAD = 4  # Bond Type
    NOESY = 6  # Through Space
    CARBONYL = 7  # Bond Type
    NONE = 9

    BOND_TYPES = [HSQC, INAD, CARBONYL]
    THROUGH_SPACE = [COSY, HMBC, NOESY]



class TwoDSignalManager:
    def __init__(self, oned_signal_manager):
        self.oned_signal_manager = oned_signal_manager
        self.number_signals = oned_signal_manager.number_signals
        self.twod_signals = []
        self.nmr_types = {"COSY" : ("H", "H", InteractionValues.COSY),
                          "HSQC" : ("C", "H", InteractionValues.HSQC),
                          "HMBC" : ("C", "H", InteractionValues.HMBC),
                          "NOESY": ("H", "H", InteractionValues.NOESY)
                          }
        self.defined_nmrs = []

    def add_nmr_signals(self, input_string, signal_type, signal_format):
        self.defined_nmrs.append(signal_type)
        if signal_type not in self.nmr_types:
            raise AttributeError("%s not a valid signal type" % signal_type)

        if signal_format == "peak":
            self.parse_peak_signals(input_string, signal_type)
        elif signal_format == "integral":
            raise NotImplemented("Haven't Processed 2D Integrals Yet...")
        else:
            raise AttributeError("%s not a valid format type" % signal_format)

    def parse_peak_signals(self, peaks_string, nmr_type):
        """
        Parses peaks from file to table
        """

        twod_peaks_table = parse_table(peaks_string, start_line=1)
        for peak in twod_peaks_table:
            signal_type_pair = self.nmr_types[nmr_type]
            signal1 = self.get_closest_1d_signal(peak[0], peak[1], peak[4], signal_type_pair[0])
            signal2 = self.get_closest_1d_signal(peak[2], peak[3], peak[5], signal_type_pair[1])
            if signal1 != signal2:
                self.twod_signals.append(TwoDSignal(signal1, signal2, nmr_type))

    def get_closest_1d_signal(self, ppm_shift, freq_shift, freq_width, signal_type):
        """
        Returns 1d signal with shift closest to input shift
        Rounds peaks to nearest valid coordinate, else Raises Exception
        """
        peak_width_scale = 2  # Fudge Factor
        freq_spectrometer = freq_shift/ppm_shift
        ppm_width = freq_width/freq_spectrometer
        dx = ppm_width * peak_width_scale
        closest_peak = None
        for signal in self.oned_signal_manager.signals:
            if abs(ppm_shift - signal.x_shift) < dx and signal_type == signal.signal_type:
                closest_peak = signal
                dx = abs(ppm_shift - signal.x_shift)

        if not closest_peak:
            for signal in self.oned_signal_manager.signals:
                print(ppm_shift, abs(ppm_shift - signal.x_shift), ppm_width)
            error_msg = "Invalid %s Peak: 2D Peak at %s, " \
                        "Does not correspond to 1D Peak" % (signal_type, ppm_shift)
            raise Exception(error_msg)

        return closest_peak

    def get_interaction_data(self):
        print("Defined NMRS:", self.defined_nmrs)
        interaction_matrix = self.get_interaction_matrix()
        type_array = self.get_type_array()
        shift_values = self.get_shift_values()
        return interaction_matrix, type_array, shift_values

    def get_shift_values(self):
        return [signal.x_shift for signal in self.oned_signal_manager.signals]

    def get_type_array(self):
        return [signal.signal_type for signal in self.oned_signal_manager.signals]

    def get_interaction_matrix(self):
        interaction_matrix = np.zeros((self.number_signals, self.number_signals), dtype=np.int)
        np.fill_diagonal(interaction_matrix, InteractionValues.NONE)
        for signal in self.twod_signals:
            for interaction in signal.get_interactions_list():
                if interaction[0] != interaction[1]:
                    interaction_matrix[interaction[0]][interaction[1]] = self.nmr_types[signal.nmr_type][2]
                    interaction_matrix[interaction[1]][interaction[0]] = self.nmr_types[signal.nmr_type][2]
        for index, row in enumerate(interaction_matrix):
            if list(row).count(InteractionValues.HSQC) > 1:
                if self.oned_signal_manager.signals[index].signal_type == "H":
                    raise Exception("H Signal with Shift %s has multiple HSQC interactions" %
                                    self.oned_signal_manager.signals[index].x_shift)
        interaction_matrix = self.infer_inadequate(interaction_matrix)
        return interaction_matrix

    def infer_inadequate(self, interaction_matrix):
        """
        Who needs a description when you can have a picture?
                         COSY
                      H-------H
                HSQC  |       |  HSQC
                      C< INAD >C
          implies  =>  INAD bond between CC
        """
        for i in range(len(interaction_matrix)):
            cosy_i = [x for x, y in enumerate(interaction_matrix[i]) if y == InteractionValues.COSY]
            hsqc_i = [x for x, y in enumerate(interaction_matrix[i]) if y == InteractionValues.HSQC]
            b_cs = [x for x in cartesian(*[cosy_i, hsqc_i])]
            for b_c in b_cs:
                b = b_c[0]
                c = b_c[1]
                d_s = [x for x, y in enumerate(interaction_matrix[b]) if y == InteractionValues.HSQC]
                b_ds = [x for x in cartesian([c], d_s)]
                for b_d in b_ds:
                    if b_d[0] != b_d[1]:
                        interaction_matrix[b_d[0]][b_d[1]] = InteractionValues.INAD
        return interaction_matrix


class TwoDSignal:
    def __init__(self, signal1, signal2, nmr_type):
        self.x_shift = round(float(signal1.x_shift), 2)
        self.y_shift = round(float(signal2.x_shift), 2)
        self.x_signal_numbers = signal1.signal_numbers
        self.y_signal_numbers = signal2.signal_numbers
        self.nmr_type = nmr_type

    def get_interactions_list(self):
        return list(cartesian(*[self.x_signal_numbers, self.y_signal_numbers]))

def plot_scatter(peaks):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.grid(True, linestyle='-', color='0.75')
    x_coords = [x[0] for x in peaks]
    y_coords = [x[1] for x in peaks]
    ax.scatter(x_coords, y_coords, s=20, marker='o')
    plt.show()
