import numpy as np

def read_signal_data(path):
    signal_file = open(path)
    signal_string = signal_file.read()
    signal_file.close()
    return signal_string

def parse_signal_string(signal_string):
    signal_string = signal_string.replace("\t"," ")
    signal_list = signal_string.splitlines()
    signal_list = [x.split(" ") for x in signal_list]
    signal_list = [[float(y) for y in x[0:3]] for x in signal_list]
    return np.array(signal_list)


def main():
    signal_string = read_signal_data("resources/hydrogen_integration_data.txt")
    parsed_signal_list = parse_signal_string(signal_string)

    

    print(parsed_signal_list)

if __name__=="__main__":
    main()
