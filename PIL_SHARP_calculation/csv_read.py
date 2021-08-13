"""
This script is written for reading all the .csv files from SHARP_with_PIL.py and consolidate them into one .csv file,
  following which we would analyze correlation between the SHARPs and flare cursting.

input: ALl .csv files including SHARP data that are computed from original data and PIL-masked data.
output: Two .csv file, one that contains all the SHARP data for one active region from original data
          while the other contains SHARP from PIL-masked data.
"""

import os, numpy as np, csv

def dict2list(dict):
    """
    Turn the dicts from .csv files to list for convenience. Normally, 'dict' is the 'row' in reader.
    :param dict: dict data directly read from .csv files
    :return:a list containing all the dict.values() in the form of float
    """
    values = dict.values()
    try:
        list = [float(x) for x in values]
    except:
        list = [x for x in values]
    return list



def main():
    path = r'res/SHARP_csv'
    list_path = os.listdir(path)

    csv_pilmasked = [x for x in list_path if 'pilmasked' in x]
    csv_pilmasked.sort()
    csv_original_cal = [x for x in list_path if 'original_calculated' in x]
    csv_original_cal.sort()

    count_file = len(csv_pilmasked)

    # Read all the .csv files and save them in a list
    sharp_name = ['totusjh', 'totpot', 'totusjz', 'absnjzh', 'savncpp', 'usflux', 'meanpot', 'r_value', \
                  'meanshr', 'meangam', 'meangbt', 'meangbz', 'meangbh', 'meanjzh', 'meanjzd', 'meanalp', 'date']

    # index change means

    SHARP_parameters = [[] for _ in range(len(sharp_name))]

    # Read all the pil_masked .csv files and put them in a variable.
    for k in range(count_file):
        with open(r'res/SHARP_csv/' + csv_pilmasked[k], newline = '') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                values_list = dict2list(row)
                for j in range(len(SHARP_parameters)):
                    SHARP_parameters[j].append(values_list[j])

    # Save the SHARP_list into one .csv file
    with open(r'res/SHARP_csv/csv_all/SHARP_pilmasked.csv', 'w', newline = '') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(sharp_name)
        for i in range(len(SHARP_parameters[0])):
            temp = [x[i] for x in SHARP_parameters]
            spamwriter.writerow(temp)

    SHARP_parameters = [[] for _ in range(len(sharp_name))]
    for k in range(count_file):
        with open (r'res/SHARP_csv/' + csv_original_cal[k], newline = '') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                values_list = dict2list(row)
                for j in range(len(SHARP_parameters)):
                    SHARP_parameters[j].append(values_list[j])

    with open(r'res/SHARP_csv/csv_all/SHARP_original_calculated.csv', 'w', newline = '') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(sharp_name)
        for i in range(len(SHARP_parameters[0])):
            temp = [x[i] for x in SHARP_parameters]
            spamwriter.writerow(temp)





    return 0

if __name__ == "__main__":
    main()
