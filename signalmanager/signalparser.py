import numpy as np


def parse_integration_string(integration_string):
    integration_string = integration_string.replace("\t", " ")
    integration_table = integration_string.splitlines()
    integration_table = [x.split(" ") for x in integration_table]
    integration_table = [[float(y) for y in x[0:3]] for x in integration_table]
    integration_table.sort()
    return np.array(integration_table)


def parse_peaks_string(twod_peaks_string):
    peak_string_list = twod_peaks_string.splitlines()
    peak_string_list = peak_string_list[1:] # Removes the Title Line
    peak_points = [[float(y) for y in x.split("\t")[0:2]] for x in peak_string_list]
    if [] in peak_points:
        peak_points.remove([])
    return peak_points