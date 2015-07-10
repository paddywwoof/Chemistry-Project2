"""
"""

def round_to_nearest(number_list, x):
    """
    Rounds x to nearest number in number_list
    """

    number_list.sort()
    for index, number in enumerate(number_list):
        if index < len(number_list)-1:
            if between(x, [number_list[index], number_list[index+1]]):
                return closest(x, [number_list[index], number_list[index+1]])
        else:
            return number_list[-1]

def between(x, interval):
    """
    Returns True if x is in interval, else returns False
    """
    if interval[0] <= x < interval[1]:
        return True
    else:
        return False

def closest(x, options):
    """
    Returns number in options x is closest to
    """
    dx = 100
    best_option=None
    for option in options:
        if abs(x-option) < dx:
            best_option = option
            dx = abs(x-option)
    return best_option

def merge_duplicates(peak_points):
    """
    """
    for peak_point1 in peak_points:
        for peak_point2 in peak_points:
            if peak_point1[0:2] == peak_point2[0:2] and peak_point1 != peak_point2:
                return merge_duplicates([x for x in peak_points if x not in [peak_point1,peak_point2]]+ [[ peak_point1[0], peak_point1[1], peak_point1[2] + peak_point2[2] ]])
    return peak_points


def remove_asymmetric(cosy_interactions):
    for interaction in cosy_interactions:
        if [interaction[1], interaction[0]] not in cosy_interactions:
            return remove_asymmetric([x for x in cosy_interactions if x not in [interaction, [interaction[1], interaction[0]] ]])
    return cosy_interactions

