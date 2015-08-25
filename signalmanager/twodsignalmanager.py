import numpy as np
import matplotlib.pyplot as plt
from filemanager import readfile
from itertools import product as cartesian
from . import OneDSignalManager
from .signalparser import parse_peaks_string

np.set_printoptions(suppress=True)

class TwoDSignalManager:
    def __init__(self, oned_signal_manager):
        self.oned_signal_manager = oned_signal_manager
        self.number_signals = oned_signal_manager.number_signals
        self.twod_signals = []
        self.nmr_types = {"COSY": ("H", "H", 1),
                          "HSQC": ("C", "H", 2),
                          "HMBC": ("C", "H", 3),
                          "NOESY": ("H","H", 6)
                          }
        self.peak_bounds = {"H": 8.5, "C": 225}
        self.shift_errors = {"H": 0.05, "C": 0.25}
        self.defined_nmrs = []

    def add_nmr_signals(self, nmr_type, twod_peaks_string):
        self.defined_nmrs.append(nmr_type)
        if nmr_type not in self.nmr_types:
            raise AttributeError("%s is not a valid NMR type" % nmr_type)
        self.parse_signals(twod_peaks_string, nmr_type)

    def parse_signals(self, twod_peaks_string, nmr_type):
        """
        Parses peaks from file to table
        Remove small peaks and those outside the spectrum
        Rounds peaks to nearest valid coordinate, or removes if too far

        If COSY, removes asymmetric peaks
        """
        peaks_table = parse_peaks_string(twod_peaks_string)
        peaks_table = self.clean_peaks_table(peaks_table, nmr_type)
        for peak in peaks_table:


            signal_type_pair = self.nmr_types[nmr_type]
            signal1 = self.get_closest_1d_signal(peak[0], self.oned_signal_manager, signal_type_pair[0])
            signal2 = self.get_closest_1d_signal(peak[1], self.oned_signal_manager, signal_type_pair[1])
            if signal1 and signal2 and signal1 != signal2:
                self.twod_signals.append(TwoDSignal(signal1, signal2, nmr_type))
            else:
                print("%s Peak Removal: %s Varies by more than 0.1 from any peak" % ( nmr_type, str(peak[0:2])) )
        if nmr_type == "COSY":
            self.twod_signals = self.remove_asymmetric(self.twod_signals)

    def get_closest_1d_signal(self, x_shift, oned_signal_manager, signal_type):
        """
        Returns 1d signal with shift closest to input shift
        """
        dx = self.shift_errors[signal_type]
        closest_peak = None
        for signal in oned_signal_manager.signals:
            if abs(x_shift - signal.x_shift) < dx and signal_type == signal.signal_type:
                closest_peak = signal
                dx = abs(x_shift - signal.x_shift)
        return closest_peak

    def remove_asymmetric(self, twod_signals):
        """
        Removes all twod signals that are not symmetric: Used for COSY spectrum
        """
        interactions_list = [[x.x_signal_numbers, x.y_signal_numbers] for x in twod_signals]
        for signal in twod_signals:
            if [signal.y_signal_numbers, signal.x_signal_numbers] not in interactions_list:
                print("COSY Peak Removal: Signal %s at point %s is asymmetric" % (str([signal.x_signal_numbers, signal. y_signal_numbers]), str([signal.x_shift, signal.y_shift])))
                twod_signals.remove(signal)
                return self.remove_asymmetric(twod_signals)
        return twod_signals

    def clean_peaks_table(self, peaks_table, nmr_type):
        """
        Removes peaks with low area or that lie outside of the nmr2 spectrum
        """
        x_signal_type = self.nmr_types[nmr_type][0]
        y_signal_type = self.nmr_types[nmr_type][1]
        for peak in peaks_table:
            """
            if peak[2] < 30 and nmr_type == "COSY":
                print("%s Peak Removal: Peak %s is too small" % (nmr_type, peak))
                peaks_table.remove(peak)
                return self.clean_peaks_table(peaks_table, nmr_type)
            """
            if not(0 < peak[0] < self.peak_bounds[x_signal_type]) or not(0 < peak[1] < self.peak_bounds[y_signal_type]):
                print("%s Peak Removal: Peak %s does not fall in correct boundary" % (nmr_type, peak))
                peaks_table.remove(peak)
                return self.clean_peaks_table(peaks_table, nmr_type)
        return peaks_table


    def get_interaction_matrix(self):
        print("Defined NMRS:", self.defined_nmrs)
        interaction_matrix = np.zeros((self.number_signals, self.number_signals), dtype=np.int)
        np.fill_diagonal(interaction_matrix, 9)
        for signal in self.twod_signals:
            for interaction in signal.get_interactions_list():
                if interaction[0] != interaction[1]:
                    interaction_matrix[interaction[0]][interaction[1]] = self.nmr_types[signal.nmr_type][2]
                    interaction_matrix[interaction[1]][interaction[0]] = self.nmr_types[signal.nmr_type][2]
        interaction_matrix = self.infer_inadequate(interaction_matrix)
        return interaction_matrix

    def infer_inadequate(self, interaction_matrix):
        """
        Who needs a description when you can have a picture?
                       cosy
                      H----H
                hsqc  |    |  hsqc
                      C<<>>C
          implies  =>  inad bond between CC
        """

        for i in range(len(interaction_matrix)):
            cosy_i = [x for x, y in enumerate(interaction_matrix[i]) if y == 1]
            hsqc_i = [x for x, y in enumerate(interaction_matrix[i]) if y == 2]
            b_cs = [x for x in cartesian(*[cosy_i, hsqc_i])]
            for b_c in b_cs:
                b = b_c[0]
                c = b_c[1]
                d_s = [x for x, y in enumerate(interaction_matrix[b]) if y == 2]
                b_ds = [x for x in cartesian([c], d_s)]
                for b_d in b_ds:
                    if b_d[0] != b_d[1]:
                        interaction_matrix[b_d[0]][b_d[1]] = 4
        return interaction_matrix

    def plot_twod_signals(self, nmr_type):
        plot_scatter([[signal.x_shift, signal.y_shift] for signal in self.twod_signals if signal.nmr_type == nmr_type])

    def get_interaction_data(self):
        type_array = []
        signals = self.oned_signal_manager.signals
        for signal in signals:
            type_array.append(signal.signal_type)
        shift_values = self.get_shift_values()
        interaction_matrix = self.get_interaction_matrix()
        return interaction_matrix, type_array, shift_values

    def get_shift_values(self):
        shift_values = []
        for signal in self.oned_signal_manager.signals:
            shift_values.append(signal.x_shift)
        return shift_values

class TwoDSignal:
    def __init__(self, signal1, signal2, nmr_type):
        self.x_shift = round(float(signal1.x_shift), 2)
        self.y_shift = round(float(signal2.x_shift), 2)
        self.x_signal_numbers = signal1.signal_numbers
        self.y_signal_numbers = signal2.signal_numbers
        self.nmr_type = nmr_type

    def get_interactions_list(self):
        return [x for x in cartesian(*[self.x_signal_numbers, self.y_signal_numbers])]

def plot_scatter(peaks):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.grid(True, linestyle='-', color='0.75')
    x_coords = [x[0] for x in peaks]
    y_coords = [x[1] for x in peaks]
    ax.scatter(x_coords, y_coords, s=20, marker='o')
    plt.show()