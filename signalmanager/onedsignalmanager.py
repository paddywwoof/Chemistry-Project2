from .signalparsers import parse_table, parse_1d_peak_list

class OneDSignalManager:
    def __init__(self):
        self.number_signals = 0
        self.signals = []

    def add_nmr_signals(self, input_string, signal_type, signal_format):

        if signal_type not in ["H", "C"]:
            raise AttributeError("%s not a valid signal type" % signal_type)
        if signal_format == "integral":
            self.parse_integration_signals(input_string, signal_type)
        elif signal_format == "peak":
            self.parse_peak_signals(input_string, signal_type)
        else:
            raise AttributeError("%s not a valid format type" % signal_format)

    def parse_peak_signals(self, peak_string, signal_type):
        peak_table = parse_1d_peak_list(peak_string, start_line=1)
        for row in peak_table:
            if row[1] == "Solvent":
                print("Solvent Peak at shift %s" % row[0])
            else:
                new_signal = OneDSignal(row[0], row[0], 1, signal_type, self.number_signals)
                self.signals.append(new_signal)
                self.number_signals += new_signal.multiplicity

    def parse_integration_signals(self, integration_string, signal_type):
        integration_table = parse_table(integration_string, start_line=0)
        for peak in integration_table:
            shift1 = peak[0]
            shift2 = peak[1]
            magnitude = peak[2]
            new_signal = OneDSignal(shift1, shift2, magnitude, signal_type, self.number_signals)
            for x in range(new_signal.multiplicity):
                self.signals.append(new_signal)
            self.number_signals += new_signal.multiplicity

    def remove_signal(self, signal):
        self.signals = [x for x in self.signals if x.x_shift != signal.x_shift]
        self.number_signals -= signal.multiplicity
        self.update_signal_numbers()

    def update_signal_numbers(self):
        start_index = 0
        for signal in self.signals:
            signal.update_signal_numbers(start_index)
            start_index += signal.multiplicity

    def add_signal(self, x_shift1, magnitude, signal_type):
        new_signal = OneDSignal(x_shift1, x_shift1, magnitude, signal_type, 0, False)
        for x in range(new_signal.multiplicity):
            self.signals.append(new_signal)
            self.number_signals += new_signal.multiplicity
        self.update_signal_numbers()


class OneDSignal:
    def __init__(self, x_shift1, x_shift2, magnitude, signal_type, start_index, round_shift=True):
        self.x_shift = 0.5*x_shift1 + 0.5*x_shift2
        if round_shift:
            self.x_shift = round(self.x_shift, 2)
        self.signal_type = signal_type
        self.multiplicity = int(round(magnitude))
        self.signal_numbers = list(range(start_index, start_index + self.multiplicity))

    def update_signal_numbers(self, start_index):
        self.signal_numbers = list(range(start_index, start_index + self.multiplicity))
