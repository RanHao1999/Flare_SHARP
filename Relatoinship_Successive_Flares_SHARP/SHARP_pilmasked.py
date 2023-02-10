"""
This code is for calculating SHARP keys manually.
input: HMI data from: http://jsoc.stanford.edu/ajax/lookdata.html
       data should be put in the 'data' folder.
       We may use the following segments:

       [example filename]                 --> [description]
        hmi.sharp_cea_*.Br.fits            --> radial component of the magnetic field vector
        hmi.sharp_cea_*.Bt.fits            --> theta-component of the magnetic field vector
        hmi.sharp_cea_*.Bp.fits            --> phi-component of the magnetic field vector
        hmi.sharp_cea_*.Br_err.fits        --> error in radial component of the magnetic field vector
        hmi.sharp_cea_*.Bt_err.fits        --> error in theta-component of the magnetic field vector
        hmi.sharp_cea_*.Bp_err.fits        --> error in phi-component of the magnetic field vector
        hmi.sharp_cea_*.conf_disambig.fits --> bits indicate confidence levels in disambiguation result
        hmi.sharp_cea_*.bitmap.fits        --> bits indicate result of automatic detection algorithm
        hmi.sharp_cea_*.magnetogram.fits   --> line-of-sight component of the magnetic field

output: .txt file as a table, including the array of all SHARP parameters needed in time series.
         16 SHARP parameters are included:
         TOTUSJH: Total unsigned current helicity  \ H_{total} = \Sigma abs(B_z \cdot J_z)

original author: Monica Bobra in 1 May 2015

modified by: Ran Hao in 09 Aug 2021

Modification made for calculating a batch of data.
"""
# import modules
import os, csv
import pandas as pd
import sunpy, sunpy.map, scipy, numpy as np, sys, math, argparse, normal_funcs as nf
import concurrent.futures
from sharpkey_pilmasked import compute_abs_flux, compute_bh, compute_gamma,\
    compute_bt, computeBtderivative, computeBhderivative, computeBzderivative, computeJz, computeJzmoments,\
    computeAlpha, computeHelicity, computeSumAbsPerPolarity, greenpot, computeFreeEnergy, computeShearAngle, computeR, computeLOSderivative, compute_abs_flux_los

radsindeg = np.pi/180.
munaught  = 0.0000012566370614
SHARP_name = ['totusjh', 'totpot', 'totusjz', 'absnjzh', 'savncpp', 'usflux', 'meanpot', 'r_value',\
              'meanshr', 'meangam', 'meangbt', 'meangbz', 'meangbh', 'meanjzh', 'meanjzd', 'meanalp', 'date']

def get_date(bz):
    time_str = bz.date.value
    return time_str

def get_data(input):
    """function: get_data

    This function reads the appropriate data and metadata.
    """
    [file_bz, file_by, file_bx, file_bz_err, file_by_err, file_bx_err, file_conf_disambig, file_bitmap, file_los, count, pilmask] = input

    try:
        bz_map = sunpy.map.Map(file_bz)
    except:
        print("Could not open the bz fits file")
        sys.exit(1)

    try:
        by_map = sunpy.map.Map(file_by)
    except:
        print("Could not open the by fits file")
        sys.exit(1)

    try:
        bx_map = sunpy.map.Map(file_bx)
    except:
        print("Could not open the bx fits file")
        sys.exit(1)

    try:
        bz_err_map = sunpy.map.Map(file_bz_err)
    except:
        print("Could not open the bz_err fits file")
        sys.exit(1)

    try:
        by_err_map = sunpy.map.Map(file_by_err)
    except:
        print("Could not open the by_err fits file")
        sys.exit(1)

    try:
        bx_err_map = sunpy.map.Map(file_bx_err)
    except:
        print("Could not open the bx_err fits file")
        sys.exit(1)

    try:
        conf_disambig_map = sunpy.map.Map(file_conf_disambig)
    except:
        print("Could not open the conf_disambig fits file")
        sys.exit(1)

    try:
        bitmap_map = sunpy.map.Map(file_bitmap)
    except:
        print("Could not open the bitmap fits file")
        sys.exit(1)

    try:
        los_map = sunpy.map.Map(file_los)
    except:
        print("Could not open the LoS fits file")
        sys.exit(1)

    #get date
    date = get_date(bz_map)

    # get array data
    bz = bz_map.data
    by = by_map.data
    bx = bx_map.data
    bz_err = bz_err_map.data
    by_err = by_err_map.data
    bx_err = bx_err_map.data
    conf_disambig = conf_disambig_map.data
    bitmap = bitmap_map.data
    los = los_map.data

    # get metadata
    header = bz_map.meta

    # get fits header key information
    rsun_ref = header['rsun_ref']
    dsun_obs = header['dsun_obs']
    rsun_obs = header['rsun_obs']
    cdelt1 = header['cdelt1']

    # Note that the value of CDELT1 in hmi.sharp_cea_720s is in units of degrees per pixel.
    # The following calculation converts CDELT1 into arcseconds.
    # Therefore the variable cdelt1_arcseconds is in units of arseconds per pixel.
    # For an explanation of this formula, see cdelt1_arcsec.pdf in this same directory.
    cdelt1_arcsec = (math.atan((rsun_ref * cdelt1 * radsindeg) / (dsun_obs))) * (1 / radsindeg) * (3600.)

    # get dimensions
    nx = bz.shape[1]
    ny = bz.shape[0]

    # Create an error array to calculate uncertainties in the keywords dervied from
    # line-of-sight data. Liu et al. (2012) [DOI: 10.1007/s11207-012-9976-x] determined
    # the median noise in the HMI full-disk magnetic field maps is 6.4 Mx cm^(−2)
    # (see Figure 2). We will assume this noise value is homogeneous throughout the disk to
    # estimate the error in the keyword quantities. Here, 1 Gauss = 1 Mx cm^(−2).
    los_err = np.ndarray(shape=(ny, nx), dtype=float)
    los_err.fill(6.4)

    # flip the sign of by
    by_flipped = -1.0 * (np.array(by))

    return [bz, by_flipped, bx, bz_err, by_err, bx_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs,
            cdelt1_arcsec, los, los_err, date, count, pilmask]

def main_run(input):
    [bz, by, bx, bz_err, by_err, bx_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, los, los_err, date, count, pilmask] = get_data(input)

    if len(pilmask) > 0:
        print("Following results SHARP parameters that're pilmasked")
        print('These are the active region summary parameters calculated from the vector magnetic field data:')
        # compute the total unsigned flux and associated errors
        mean_vf, mean_vf_err, count_mask = compute_abs_flux(bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, pilmask)
        print('USFLUX ', mean_vf, 'Mx')
        print('ERRVF', mean_vf_err, 'Mx')
        print('CMASK', count_mask, 'pixels')

        # compute the horizontal component of the magnetic field and associated errors
        horiz = compute_bh(bx, by, bz, bx_err, by_err, bz_err, conf_disambig, bitmap, nx, ny, pilmask)
        bh, bh_err = horiz[0], horiz[1]

        # compute the shear angle and associated errors
        mean_gamma, mean_gamma_err = compute_gamma(bx, by, bz, bh, bz_err, bh_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, pilmask)
        print('MEANGAM ', mean_gamma, 'degree')
        print('ERRGAM ', mean_gamma_err, 'degree')

        # compute the total magnetic field vector and associated errors
        total = compute_bt(bx, by, bz, bx_err, by_err, bz_err, conf_disambig, bitmap, nx, ny, pilmask)
        bt, bt_err = total[0], total[1]

        # compute the field gradients and associated errors
        mean_derivative_bt, mean_derivative_bt_err = computeBtderivative(bt, bt_err, nx, ny, conf_disambig, bitmap, pilmask)
        print('MEANGBT ', mean_derivative_bt, 'G * Mm^(-1)')
        print('ERRBT ', mean_derivative_bt_err, 'G * Mm^(-1)')

        mean_derivative_bh, mean_derivative_bh_err = computeBhderivative(bh, bh_err, nx, ny, conf_disambig, bitmap, pilmask)
        print('MEANGBH ', mean_derivative_bh, 'G * Mm^(-1)')
        print('ERRBH ', mean_derivative_bh_err, 'G * Mm^(-1)')

        mean_derivative_bz, mean_derivative_bz_err = computeBzderivative(bz, bz_err, nx, ny, conf_disambig, bitmap, pilmask)
        print('MEANGBZ ', mean_derivative_bz, 'G * Mm^(-1)')
        print('ERRBZ ', mean_derivative_bz_err, 'G * Mm^(-1)')

        # compute the vertical current and associated errors
        current = computeJz(bx, by, bx_err, by_err, conf_disambig, bitmap, nx, ny, pilmask)
        jz, jz_err, derx, dery = current[0], current[1], current[2], current[3]

        # compute the moments of the vertical current density and associated errors
        mean_jz, mean_jz_err, us_i, us_i_err = computeJzmoments(jz, jz_err, derx, dery, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, munaught, pilmask)
        print('MEANJZD ', mean_jz, 'mA * m^(−2)')
        print('ERRJZ ', mean_jz_err, 'mA * m^(−2)')
        print('TOTUSJZ ', us_i, 'A')
        print('ERRUSI', us_i_err, 'A')

        # compute the twist parameter, alpha, and associated errors
        mean_alpha, mean_alpha_err = computeAlpha(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, pilmask)
        print('MEANALP ', mean_alpha, 'Mm^(-1)')
        print('ERRALP ', mean_alpha_err, 'Mm^(-1)')

        # compute the moments of the current helicity and associated errors
        mean_ih, mean_ih_err, total_us_ih, total_us_ih_err, total_abs_ih, total_abs_ih_err = computeHelicity(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, pilmask)
        print('MEANJZH ', mean_ih, 'G2 * m^(−1)')
        print('ERRMIH ', mean_ih_err, 'G2 * m^(−1)')
        print('TOTUSJH ', total_us_ih, 'G2 * m^(−1)')
        print('ERRTUI ', total_us_ih_err, 'G2 * m^(−1)')
        print('ABSNJZH ', total_abs_ih, 'G2 * m^(−1)')
        print('ERRTAI ', total_abs_ih_err, 'G2 * m^(−1)')

        # compute the sum of the absolute value per polarity and associated errors
        totaljz, totaljz_err = computeSumAbsPerPolarity(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, munaught, pilmask)
        print('SAVNCPP ', totaljz, 'A')
        print('ERRJHT ', totaljz_err, 'A')

        # compute a numerical model of the potential field (it has no errors, as the theoretical values are exact)
        potential = greenpot(bz, nx, ny)
        bpx, bpy = potential[0], potential[1]

        # compute the energy stored in the magnetic field and its associated errors
        meanpot, meanpot_err, totpot, totpot_err = computeFreeEnergy(bx_err, by_err, bx, by, bpx, bpy, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, conf_disambig, bitmap, pilmask)
        print('MEANPOT ', meanpot, 'erg * cm^(−3)')
        print('ERRMPOT ', meanpot_err, 'erg * cm^(−3)')
        print('TOTPOT ', totpot, 'erg * cm^(−1)')
        print('ERRTPOT ', totpot_err, 'erg * cm^(−1)')

        # compute the degree to which the observed field is sheared and its associated errors
        meanshear_angle, meanshear_angle_err, area_w_shear_gt_45 = computeShearAngle(bx_err, by_err, bz_err, bx, by, bz, bpx, bpy, nx, ny, conf_disambig, bitmap, pilmask)
        print('MEANSHR ', meanshear_angle, 'degree')
        print('ERRMSHA ', meanshear_angle_err, 'degree')
        print('SHRGT45 ', area_w_shear_gt_45, 'as a percentage')

        print('===============================')
        print('These are the active region summary parameters calculated from the line-of-sight magnetic field data:')
        # compute the gradient-weighted neutral line length
        Rparam, Rparam_err = computeR(los, los_err, nx, ny, cdelt1_arcsec, pilmask)
        print('R_VALUE ', Rparam, 'Mx')
        print('The error in R_VALUE is', Rparam_err)

        # compute mean gradient of the line-of-sight field
        mean_derivative_blos, mean_derivative_blos_err = computeLOSderivative(los, los_err, nx, ny, bitmap, rsun_ref, rsun_obs, cdelt1_arcsec, pilmask)
        print('MEANGBL ', mean_derivative_blos, 'G * Mm^(-1)')
        print('The error in MEANGBL is', mean_derivative_blos_err)

        # compute the total unsigned flux using the line of sight field
        mean_vf, mean_vf_err, count_mask = compute_abs_flux_los(los, los_err, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, pilmask)
        print('USFLUXL ', mean_vf, 'Mx')
        print('The error in USFLUX is', mean_vf_err)
        print('CMASKL', count_mask, 'pixels')

        print('Note that the calculation for R_VALUE uses a slightly different method than applied for the hmi.sharp*_720s series. The results, however, should be identical or within a log(R) value of 0.1. ')
        print('All the other keyword calculations use an identical method, and the results are identical. ')

        SHARP_row = [total_us_ih, totpot, us_i, total_abs_ih, totaljz, mean_vf, meanpot, Rparam, meanshear_angle, mean_gamma, mean_derivative_bt, mean_derivative_bz, mean_derivative_bh, mean_ih, mean_jz, mean_alpha, date]
        with open('res/SHARP_csv/SHARP_pilmasked_' + str('%04d' % count) + '.csv', 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile)
            spamwriter.writerow(SHARP_name)
            spamwriter.writerow(SHARP_row)

    else:
        print('empty PIL')
        SHARP_row = [0 for _ in range(len(SHARP_name))]
        SHARP_row[-1] = date
        with open('res/SHARP_csv/SHARP_pilmasked_' + str('%04d' % count) + '.csv', 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile)
            spamwriter.writerow(SHARP_name)
            spamwriter.writerow(SHARP_row)

    return 0

def main():
    # Read all the .fits data and put them in a
    data_dir = 'data/JSOC'
    file_lst = os.listdir(data_dir)

    pil_dir = 'res/coord_pil'
    file_pils = os.listdir(pil_dir)

    bz_lst = [x for x in file_lst if 'Br' in x and 'Br_err' not in x]
    by_lst = [x for x in file_lst if 'Bt' in x and 'Bt_err' not in x]
    bx_lst = [x for x in file_lst if 'Bp' in x and 'Bp_err' not in x]
    bz_err_lst = [x for x in file_lst if 'Br_err' in x]
    by_err_lst = [x for x in file_lst if 'Bt_err' in x]
    bx_err_lst = [x for x in file_lst if 'Bp_err' in x]
    mask_lst = [x for x in file_lst if 'conf_disambig' in x]
    bitmask_lst = [x for x in file_lst if 'bitmap' in x]
    los_lst = [x for x in file_lst if 'magnetogram' in x]
    pil_lst = [x for x in file_pils if '.txt' in x]

    bz_lst.sort()
    by_lst.sort()
    bx_lst.sort()
    bz_err_lst.sort()
    by_err_lst.sort()
    bx_err_lst.sort()
    mask_lst.sort()
    bitmask_lst.sort()
    los_lst.sort()
    pil_lst.sort()
    count_file = int(len(bz_lst))
    input = []

    for count in range(count_file):
        file_bz = r'data/JSOC/' + bz_lst[count]
        file_by = r'data/JSOC/' + by_lst[count]
        file_bx = r'data/JSOC/' + bx_lst[count]
        file_bz_err = r'data/JSOC/' + bz_err_lst[count]
        file_by_err = r'data/JSOC/' + by_err_lst[count]
        file_bx_err = r'data/JSOC/' + bx_err_lst[count]
        file_conf_disambig = r'data/JSOC/' + mask_lst[count]
        file_bitmap = r'data/JSOC/' + bitmask_lst[count]
        file_los = r'data/JSOC/' + los_lst[count]
        file_pil = r'res/coord_pil/' + pil_lst[count]

        try:
            pil_table = pd.read_table(file_pil)
            pilmask = nf.dataframe2matrix(pil_table)
            pilmask = pilmask.tolist()
        except:
            pilmask = []

        res = [file_bz, file_by, file_bx, file_bz_err, file_by_err, file_bx_err, file_conf_disambig, file_bitmap, file_los]
        res.append(count)
        res.append(pilmask)
        input.append(res)

    with concurrent.futures.ProcessPoolExecutor(max_workers = 20) as pool:
        pool.map(main_run, input)


    return 0


if __name__ == "__main__":
    main()
