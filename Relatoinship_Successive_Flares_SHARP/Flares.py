"""
# Read flare data from flare list and analyze the correlation between flare bursting with original SHARP.
"""

import csv, numpy as np, datetime, math, matplotlib.pyplot as plt

import pandas as pd

from normal_funcs import dict2list, list2date0, list2date1,SHARP2float, normalization

AR_num = '11158'
sharp_name = ['totusjh', 'totpot', 'totusjz', 'absnjzh', 'savncpp', 'usflux', 'meanpot', 'r_value',\
              'meanshr', 'meangam', 'meangbt', 'meangbz', 'meangbh', 'meanjzh', 'meanjzd', 'meanalp', 'date']

def flare_level2index(flare_level):
    # C1.0 = 0.0, M1.0 = 1.0, X1.0 = 2
    index = 0

    if 'C' in flare_level:
        index = math.log(float(flare_level[1:]), 10)
    if 'M' in flare_level:
        index = math.log(10.0*float(flare_level[1:]), 10)
    if 'X' in flare_level:
        index = math.log(100.0*float(flare_level[1:]), 10)

    return index

def Read_flare(AR_num):
    path_flare = r'data/flare/xrt_flarecat.csv'

    with open(path_flare, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        list_row = [row for row in reader]

    list_class = [x['class'] for x in list_row if x['region'] == AR_num ]
    list_time = [x['peak'] for x in list_row if x['region'] == AR_num]

    return list_class, list_time

def read_SHARP_csv(file_csv):
    # Read .csv file that contains the SHARP data and return a list of all SHARP parameters.
    SHARP_parameters = [[] for _ in range(len(sharp_name))]
    with open(file_csv, newline = '') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            values_list = dict2list(row)
            for j in range(len(SHARP_parameters)):
                SHARP_parameters[j].append(values_list[j])
    return SHARP_parameters

def flare_envelope(flare_array):
    # Input an array!
    # Return the envelope curve of the flare_bursting.
    count = 0
    for i in range(len(flare_array)):
        if flare_array[i] != 0:
            k = (flare_array[i] - flare_array[count])/ (i - count + 1)
            for j in range(count+1, i):
                flare_array[j] = flare_array[count] + k * (j - count)
            count = i
    if count != len(flare_array):
        k = flare_array[count] / (len(flare_array) - count)
        for i in range(count+1, len(flare_array)):
            flare_array[i] = flare_array[count] - k*(i - count)
    return flare_array

def consolidate_SHARP_flare(SHARP_time, flare_list, flare_time):
    """
    Consolidate SHARP data and flare data together. Expand the flare data to get its envelope curve.
    :return: A flare list that's been expanded to the size of SHARP.
    """
    # Step1. Delete the flares that not happen in the SHARP_time range.
    time_tobe_removed = []
    class_tobe_removed = []
    for i in range(len(flare_time)):
        if flare_time[i] < SHARP_time[0] or flare_time[i] > SHARP_time[-1]:
            time_tobe_removed.append(flare_time[i])
            class_tobe_removed.append(flare_list[i])
    for i in range(len(time_tobe_removed)):
        flare_time.remove(time_tobe_removed[i])
        flare_list.remove(class_tobe_removed[i])

    # Step2. Expand the flare data to the range of SHARP using envelope curve method.
    New_flare_arr = np.zeros(len(SHARP_time))
    for i in range(len(flare_time)):
        for j in range(len(SHARP_time)):
            if flare_time[i] >= SHARP_time[j-1] and flare_time[i] <= SHARP_time[j]:
                New_flare_arr[j] = flare_list[i]
    New_flare_lst = flare_envelope(New_flare_arr)
    return New_flare_lst

def plot_together(data1, data2, time, name, note):
    # data1: SHARP, data2: Goes
    data1_s = pd.Series(data1)
    data2_s = pd.Series(data2)

    corre = round(data1_s.corr(data2_s), 3)
    fig, ax = plt.subplots()
    ax.plot(time, data1, 'black', label = name)
    ax.set_xlabel('Time')
    ax.set_ylabel(name)
    ax.set_title('Correlation = '+ str(corre))
    ax2 = ax.twinx()
    ax2.plot(time, data2, 'green', label = 'flare index')
    ax2.set_ylabel('flare index')
#    plt.show()
    fig.autofmt_xdate()
    fig.savefig(r'res/SHARP_Goes_pic/' + note + '/'+name+'_Goes.jpg')
    return corre

def write_csv_SF(SHARP, flare_env, note, Sharp_name, time):
    """
    This function write .csvfile in csv_SF which contains the SHARP data and the flare data altogether
    """
    SHARP = np.array(SHARP)
    SHARP = np.transpose(SHARP)
    SHARP.tolist()

    with open('res/SHARP_csv/csv_SF/SF_'+note+'.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(sharp_name + ['flare index', 'Time'])
        for i in range(len(SHARP)):
            spamwriter.writerow(SHARP[i].tolist() + [flare_env[i], time[i]])

    return 0

def smooth(y, box_pts):
    # smooth the pilmasked data
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def main():
    # Read flare bursting time
    listf_class, listf_time = Read_flare(AR_num)
    for i in range(len(listf_class)):
        listf_class[i] = flare_level2index(listf_class[i])

    #Read original SHARP data and time
    file_original = r'res/SHARP_csv/csv_all/SHARP_original.csv'
    SHARP_original = read_SHARP_csv(file_original)
    Time_SHARP = SHARP_original[-1]
    SHARP_original.remove(SHARP_original[-1])
    SHARP_original = SHARP2float(SHARP_original)

    #Read original-calculated data and time
    file_original_c = r'res/SHARP_csv/csv_all/SHARP_original_calculated.csv'
    SHARP_original_c = read_SHARP_csv(file_original_c)
    Time_SHARP_c = SHARP_original_c[-1]
    SHARP_original_c.remove(SHARP_original_c[-1])
    SHARP_original_c = SHARP2float(SHARP_original_c)

    # Read masked SHARP data
    file_masked = r'res/SHARP_csv/csv_all/SHARP_pilmasked.csv'
    SHARP_masked = read_SHARP_csv(file_masked)
    Time_SHARP_masked = SHARP_masked[-1]
    SHARP_masked.remove(SHARP_masked[-1])
    SHARP_masked = SHARP2float(SHARP_masked)

    # SHARP_masked smoothed here
    SHARP_masked_smoothed = [smooth(x, 5) for x in SHARP_masked]

    # Turn time variables to datetime object
    listf_time = list2date1(listf_time)
    Time_SHARP = list2date0(Time_SHARP)
    Time_SHARP_c = list2date0(Time_SHARP_c)
    Time_SHARP_masked = list2date0(Time_SHARP_masked)

    # Get the envelope curve of flares.
    flare_env = consolidate_SHARP_flare(Time_SHARP, listf_class, listf_time)  # flare_env is the envelop curve of flare bursting and it's the size of SHARP.
    flare_env_c = consolidate_SHARP_flare(Time_SHARP_c, listf_class, listf_time)
    flare_env_masked = consolidate_SHARP_flare(Time_SHARP_masked, listf_class, listf_time)

    # Calculate the correlation between original SHARP parameters and flare envelope curve.
    correlation_lst = []
    flare_series = pd.Series(flare_env)
    for i in range(len(SHARP_original)):
        So_series = pd.Series(np.array(SHARP_original[i]))
        Sm_series = pd.Series(np.array(SHARP_masked[i]))
        Soc_series = pd.Series(np.array(SHARP_original_c[i]))
        Sms_series = pd.Series(np.array(SHARP_masked_smoothed[i]))


        correlation_lst.append([sharp_name[i], So_series.corr(flare_series), Sm_series.corr(flare_series), Soc_series.corr(flare_series), Sms_series.corr(flare_series)])

    with open(r'res/SHARP_csv/csv_correlation/correlation1.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(['SHARP', 'original', 'pil-masked', 'original-calculated', 'pilmasked_smoothed'])
        for i in range(len(correlation_lst)):
            spamwriter.writerow(correlation_lst[i])

    for i in range(len(SHARP_original)):
        plot_together(SHARP_original[i], flare_env, Time_SHARP, sharp_name[i], 'original')

    for i in range(len(SHARP_original_c)):
        plot_together(SHARP_original_c[i], flare_env_c, Time_SHARP_c, sharp_name[i], 'original_calculated')

    for i in range(len(SHARP_masked)):
        plot_together(SHARP_masked_smoothed[i], flare_env_masked, Time_SHARP_masked, sharp_name[i], 'pilmasked_smoothed')


    sharp_name.remove(sharp_name[-1])
    write_csv_SF(SHARP_original, flare_env, 'original', sharp_name, Time_SHARP)
    write_csv_SF(SHARP_original_c, flare_env_c, 'original_calculated', sharp_name, Time_SHARP_c)
    write_csv_SF(SHARP_masked, flare_env_masked, 'pilmasked', sharp_name, Time_SHARP_masked)

    print('break here')
    return 0

if __name__ == "__main__":
    main()
































