import numpy as np
from filemanager import readfile
from .signalparser import parse_integration_string

class OneDSignalManager:
    def __init__(self):
        self.number_signals = 0
        self.signals = []

    def add_nmr_signals(self, integration_string, signal_type):
        if signal_type not in ["H", "C"]:
            raise AttributeError("%s not a valid signal type" % signal_type)
        self.parse_signals(integration_string, signal_type)

    def parse_signals(self, integration_string, signal_type):
        integration_table = parse_integration_string(integration_string)
        for peak in integration_table:
            shift1 = peak[0]
            shift2 = peak[1]
            if signal_type == "H":
                magnitude = peak[2]
            elif signal_type == "C":
                magnitude = 1
            else:
                magnitude = 0
            new_signal = OneDSignal(shift1, shift2, magnitude, signal_type, self.number_signals)
            if not(signal_type == "C" and 76.90 <= new_signal.x_shift <= 77.60):
                for x in range(new_signal.multiplicity):
                    self.signals.append(new_signal)
                self.number_signals += new_signal.multiplicity
            else:
                print("Solvent Signal Detected at %s" % new_signal.x_shift)

    def print_signals(self):
        for signal in self.signals:
            print(signal.x_shift, " "*(6-len(str(signal.x_shift))), signal.signal_type, signal.multiplicity, signal.signal_numbers)


class OneDSignal:
    def __init__(self, x_shift1, x_shift2, magnitude, signal_type, start_index):
        self.x_shift = self.get_shift(x_shift1, x_shift2)
        self.signal_type = signal_type
        self.multiplicity = int(round(magnitude))
        self.signal_numbers = self.get_signal_numbers(start_index)

    def get_signal_numbers(self, start_index):
        return [x for x in range(start_index, start_index + self.multiplicity)]

    def get_shift(self, x_shift1, x_shift2):
        x_shift_average = 0.5*x_shift1 + 0.5*x_shift2
        x_shift_rounded = round(x_shift_average, 2)
        return x_shift_rounded

def main():
    signal_manager = OneDSignalManager()
    signal_manager.add_nmr_signals('resources/oned/hydrogen_integration_data.txt', "H")
    signal_manager.add_nmr_signals('resources/oned/carbon_integration_data.txt', "C")
    signal_manager.print_signals()


if __name__ == "__main__":
    main()
