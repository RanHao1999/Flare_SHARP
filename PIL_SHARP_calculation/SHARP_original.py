"""
Read .fits files and return the SHARP parameters and write them in a .csv file.
"""

import csv, sunpy, sunpy.map, os


def main():
    path = r'data/JSOC'
    file_lst = os.listdir(path)
    file_lst.sort()
    bz_lst = [x for x in file_lst if 'Bp' in x and 'Bp_err' not in x]

    sharp_name = ['totusjh', 'totpot', 'totusjz', 'absnjzh', 'savncpp', 'usflux', 'meanpot', 'r_value',\
              'meanshr', 'meangam', 'meangbt', 'meangbz', 'meangbh', 'meanjzh', 'meanjzd', 'meanalp', 'date']

    sharp_original = [[x] for x in sharp_name]

    for count in range(len(bz_lst)):
        file_bz = r'data/JSOC/' + bz_lst[count]
        bz = sunpy.map.Map(file_bz)
        for i in range(len(sharp_original)-1):
            sharp_original[i].append(bz.meta[sharp_original[i][0]])
        sharp_original[-1].append(bz.date.value)
    with open(r'res/SHARP_csv/csv_all/SHARP_original.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile)
        for i in range(len(sharp_original[0])):
            temp = [x[i] for x in sharp_original]
            spamwriter.writerow(temp)

    return 0
if __name__ == "__main__":
    main()







