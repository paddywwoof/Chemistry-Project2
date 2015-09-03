import numpy as np


"""
def parse_integration_string(integration_string):
    integration_string = integration_string.replace("\t", " ")
    integration_table = integration_string.splitlines()
    integration_table = [x.split(" ") for x in integration_table]
    integration_table = [[float(y) for y in x[0:3]] for x in integration_table]
    integration_table.sort()
    return np.array(integration_table)
"""


def parse_peaks_string(twod_peaks_string):
    peak_string_list = twod_peaks_string.splitlines()
    peak_string_list = peak_string_list[1:] # Removes the Title Line
    peak_points = [[float(y) for y in x.split("\t")[0:2]] for x in peak_string_list]
    if [] in peak_points:
        peak_points.remove([])
    return peak_points


def parse_table(table_string, start_line=0):
    d = table_string.splitlines()[start_line:]
    e = [np.fromstring(x, sep="\t") for x in d]
    f = [x for x in e if len(x) == len(e[0])]
    return np.array(f)


def parse_1d_peak_list(list_string, start_line=0):
    d = list_string.splitlines()[start_line:]
    e = [x for x in d if len(x) > 0]
    peak_list = []
    for line in e:
        if len(line) > 0:
            line = line.replace("\t", " ")
            while "  " in line:
                line = line.replace("  ", " ")
            line_list = line.split(" ")
            peak_shift = float(line_list[0])
            peak_type = line_list[1]
            peak_list.append([peak_shift, peak_type])
    return peak_list


if __name__ == "__main__":
    """
    proton_string = open("C:\\Users\\Martin\\Desktop\\NMR Samples\\edmpc_new_format\\proton.txt").read()
    proton_table = parse_table(proton_string, start_line=0)
    print(proton_table)

    carbon_string = open("C:\\Users\\Martin\\Desktop\\NMR Samples\\edmpc_new_format\\carbon.txt").read()
    carbon_table = parse_1d_peak_list(carbon_string, start_line=1)
    [print(x) for x in carbon_table]

    cosy_string = open("C:\\Users\\Martin\\Desktop\\NMR Samples\\edmpc_new_format\\cosy.txt").read()
    cosy_table = parse_table(cosy_string, start_line=1)
    print(cosy_table)
    """
    noesy_string = open("C:\\Users\\Martin\\Desktop\\noesy.txt").read()
    noesy_table = parse_table(noesy_string, start_line=1)
    print(noesy_table)

