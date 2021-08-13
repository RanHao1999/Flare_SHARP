"""
This is the code for calculating SHARP keys with the codes provided by Monica Bobra (github @Monica Bobra).
This script contains all the functions and the functions only for calculating the SHARP keys.
"""

import sunpy.map
import scipy.ndimage
import numpy as np
import sys
import math
import argparse
from skimage.measure import block_reduce

radsindeg = np.pi/180.
munaught  = 0.0000012566370614

def compute_abs_flux(bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec):
    """function: compute_abs_flux

    This function computes the total unsigned flux in units of G/cm^2.
    It also returns the number of pixels used in this calculation in the keyword CMASK.

    To compute the unsigned flux, we simply calculate
       flux = surface integral [(vector Bz) dot (normal vector)],
            = surface integral [(magnitude Bz)*(magnitude normal)*(cos theta)].

    However, since the field is radial, we will assume cos theta = 1.
    Therefore, the pixels only need to be corrected for the projection.

    To convert G to G*cm^2, simply multiply by the number of square centimeters per pixel:
       (Gauss/pix^2)(CDELT1)^2(RSUN_REF/RSUN_OBS)^2(100.cm/m)^2
       =Gauss*cm^2
    """

    count_mask = 0
    sum = 0.0
    err = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(bz[j, i]):
                continue
            sum += abs(bz[j, i])
            err += bz_err[j, i] * bz_err[j, i]
            count_mask += 1

    mean_vf = sum * cdelt1_arcsec * cdelt1_arcsec * (rsun_ref / rsun_obs) * (rsun_ref / rsun_obs) * 100.0 * 100.0
    mean_vf_err = (np.sqrt(err)) * abs(
        cdelt1_arcsec * cdelt1_arcsec * (rsun_ref / rsun_obs) * (rsun_ref / rsun_obs) * 100.0 * 100.0)

    return [mean_vf, mean_vf_err, count_mask]


# ===========================================

def compute_bh(bx, by, bz, bx_err, by_err, bz_err, conf_disambig, bitmap, nx, ny):
    """function: compute_bh

    This function calculates B_h, the horizontal field, in units of Gauss.
    (The magnetic field has native units of Gauss since the filling factor = 1).
    """

    bh = np.zeros([ny, nx])
    bh_err = np.zeros([ny, nx])

    for j in range(ny):
        for i in range(nx):
            if (np.isnan(bx[j, i]) or np.isnan(by[j, i])):
                bh[j, i] = np.nan
                bh_err[j, i] = np.nan
                continue
            bh[j, i] = np.sqrt(bx[j, i] * bx[j, i] + by[j, i] * by[j, i])
            bh_err[j, i] = np.sqrt(
                bx[j, i] * bx[j, i] * bx_err[j, i] * bx_err[j, i] + by[j, i] * by[j, i] * by_err[j, i] * by_err[j, i]) / \
                           bh[j, i]

    return [bh, bh_err]


# ===========================================

def compute_gamma(bx, by, bz, bh, bz_err, bh_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec):
    """function: compute_gamma

    This function computes the inclination of the horizontal field (relative to the radial field).

    Error analysis calculations are done in radians (since derivatives are only true in units of radians),
    and multiplied by (180./PI) at the end for consistency in units.
    """

    count_mask = 0
    sum = 0.0
    err = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if (np.isnan(bz[j, i]) or np.isnan(bz_err[j, i]) or np.isnan(bh[j, i]) or np.isnan(bh_err[j, i]) or bz[
                j, i] == 0):
                continue
            if (bh[j, i] < 100):
                continue
            sum += abs(math.atan(bh[j, i] / abs(bz[j, i]))) * (180. / np.pi)
            err += (1 / (1 + ((bh[j, i] * bh[j, i]) / (bz[j, i] * bz[j, i])))) * (
                        1 / (1 + ((bh[j, i] * bh[j, i]) / (bz[j, i] * bz[j, i])))) * (
                               ((bh_err[j, i] * bh_err[j, i]) / (bz[j, i] * bz[j, i])) + (
                                   (bh[j, i] * bh[j, i] * bz_err[j, i] * bz_err[j, i]) / (
                                       bz[j, i] * bz[j, i] * bz[j, i] * bz[j, i])))
            count_mask += 1

    mean_gamma = sum / count_mask
    mean_gamma_err = (np.sqrt(err) / (count_mask)) * (180. / np.pi)

    return [mean_gamma, mean_gamma_err]


# ===========================================

def compute_bt(bx, by, bz, bx_err, by_err, bz_err, conf_disambig, bitmap, nx, ny):
    """function: compute_bt

    This function calculates B_t, the total field, in units of Gauss.
    (The magnetic field has native units of Gauss since the filling factor = 1).
    """

    bt = np.zeros([ny, nx])
    bt_err = np.zeros([ny, nx])

    for j in range(ny):
        for i in range(nx):
            if (np.isnan(bx[j, i]) or np.isnan(by[j, i]) or np.isnan(bz[j, i])):
                bt[j, i] = np.nan
                bt_err[j, i] = np.nan
                continue
            bt[j, i] = np.sqrt(bx[j, i] * bx[j, i] + by[j, i] * by[j, i] + bz[j, i] * bz[j, i])
            bt_err[j, i] = np.sqrt(
                bx[j, i] * bx[j, i] * bx_err[j, i] * bx_err[j, i] + by[j, i] * by[j, i] * by_err[j, i] * by_err[j, i] +
                bz[j, i] * bz[j, i] * bz_err[j, i] * bz_err[j, i]) / bt[j, i]

    return [bt, bt_err]


# ===========================================

def computeBtderivative(bt, bt_err, nx, ny, conf_disambig, bitmap):
    """function: computeBtderivative

    This function computes the derivative of the total field, or sqrt[(dB_total/dx)^2 + (dB_total/dy)^2].
    The native units in the series hmi.sharp_720s and hmi.sharp_cea_720s are in Gauss/pixel.

    Here are the steps to convert from Gauss/pixel to Gauss/Mm:
    The units of the magnetic field, or dB_total, are in Gauss.
    The units of length, i.e. dx or dy, are in pixels.
    The units of dB_total/dx or dB_total/dy = (Gauss/pix)(pix/arcsec)(arsec/meter)(meter/Mm), or
                                            = (Gauss/pix)(1/cdelt1_arcsec)(RSUN_OBS/RSUN_REF)(1000000)
                                            = Gauss/Mm
    In other words, multiply MEANGBT by a factor of (1/cdelt1_arcsec)*(RSUN_OBS/RSUN_REF)*(1000000).

    Note that cdelt1_arcsec is defined in the get_data function above.

    """

    count_mask = 0
    sum = 0.0
    err = 0.0

    derx_bt = np.zeros([ny, nx])
    dery_bt = np.zeros([ny, nx])
    err_term1 = np.zeros([ny, nx])
    err_term2 = np.zeros([ny, nx])

    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1, nx - 1):
        for j in range(0, ny):
            derx_bt[j, i] = (bt[j, i + 1] - bt[j, i - 1]) * 0.5
            err_term1[j, i] = (((bt[j, i + 1] - bt[j, i - 1]) * (bt[j, i + 1] - bt[j, i - 1])) * (
                        bt_err[j, i + 1] * bt_err[j, i + 1] + bt_err[j, i - 1] * bt_err[j, i - 1]))

    # brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0, nx):
        for j in range(1, ny - 1):
            dery_bt[j, i] = (bt[j + 1, i] - bt[j - 1, i]) * 0.5
            err_term2[j, i] = (((bt[j + 1, i] - bt[j - 1, i]) * (bt[j + 1, i] - bt[j - 1, i])) * (
                        bt_err[j + 1, i] * bt_err[j + 1, i] + bt_err[j - 1, i] * bt_err[j - 1, i]))

    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero.
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the conf_disambig and bitmap arrays.

    i = 0
    for j in range(ny):
        derx_bt[j, i] = ((-3 * bt[j, i]) + (4 * bt[j, i + 1]) - (bt[j, i + 2])) * 0.5

    i = nx - 1
    for j in range(ny):
        derx_bt[j, i] = ((3 * bt[j, i]) + (-4 * bt[j, i - 1]) - (-bt[j, i - 2])) * 0.5

    j = 0
    for i in range(nx):
        dery_bt[j, i] = ((-3 * bt[j, i]) + (4 * bt[j + 1, i]) - (bt[(j + 2), i])) * 0.5

    j = ny - 1
    for i in range(nx):
        dery_bt[j, i] = ((3 * bt[j, i]) + (-4 * bt[j - 1, i]) - (-bt[j - 2, i])) * 0.5

    # Calculate the sum only
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if ((derx_bt[j, i] + dery_bt[j, i]) == 0):
                continue
            if np.isnan(bt[j, i]):
                continue
            if np.isnan(bt[j + 1, i]):
                continue
            if np.isnan(bt[j - 1, i]):
                continue
            if np.isnan(bt[j, i - 1]):
                continue
            if np.isnan(bt[j, i + 1]):
                continue
            if np.isnan(bt_err[j, i]):
                continue
            if np.isnan(derx_bt[j, i]):
                continue
            if np.isnan(dery_bt[j, i]):
                continue
            sum += np.sqrt(derx_bt[j, i] * derx_bt[j, i] + dery_bt[j, i] * dery_bt[j, i])
            err += err_term2[j, i] / (16.0 * (derx_bt[j, i] * derx_bt[j, i] + dery_bt[j, i] * dery_bt[j, i])) + \
                   err_term1[j, i] / (16.0 * (derx_bt[j, i] * derx_bt[j, i] + dery_bt[j, i] * dery_bt[j, i]))
            count_mask += 1

    mean_derivative_bt = (sum) / (count_mask)
    mean_derivative_bt_err = (np.sqrt(err)) / (count_mask)

    return [mean_derivative_bt, mean_derivative_bt_err]


# ===========================================

def computeBhderivative(bh, bh_err, nx, ny, conf_disambig, bitmap):
    """function: computeBhderivative

    This function computes the derivative of the horizontal field, or sqrt[(dB_h/dx)^2 + (dB_h/dy)^2].
    The native units in the series hmi.sharp_720s and hmi.sharp_cea_720s are in Gauss/pixel.

    Here are the steps to convert from Gauss/pixel to Gauss/Mm:
    The units of the magnetic field, or dB_h, are in Gauss.
    The units of length, i.e. dx or dy, are in pixels.
    The units of dB_h/dx or dB_h/dy = (Gauss/pix)(pix/arcsec)(arsec/meter)(meter/Mm), or
                                    = (Gauss/pix)(1/cdelt1_arcsec)(RSUN_OBS/RSUN_REF)(1000000)
                                    = Gauss/Mm
    In other words, multiply MEANGBH by a factor of (1/cdelt1_arcsec)*(RSUN_OBS/RSUN_REF)*(1000000).

    Note that cdelt1_arcsec is defined in the get_data function above.

    """

    count_mask = 0
    sum = 0.0
    err = 0.0

    derx_bh = np.zeros([ny, nx])
    dery_bh = np.zeros([ny, nx])
    err_term1 = np.zeros([ny, nx])
    err_term2 = np.zeros([ny, nx])

    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1, nx - 1):
        for j in range(0, ny):
            derx_bh[j, i] = (bh[j, i + 1] - bh[j, i - 1]) * 0.5
            err_term1[j, i] = (((bh[j, i + 1] - bh[j, i - 1]) * (bh[j, i + 1] - bh[j, i - 1])) * (
                        bh_err[j, i + 1] * bh_err[j, i + 1] + bh_err[j, i - 1] * bh_err[j, i - 1]))

    # brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0, nx):
        for j in range(1, ny - 1):
            dery_bh[j, i] = (bh[j + 1, i] - bh[j - 1, i]) * 0.5
            err_term2[j, i] = (((bh[j + 1, i] - bh[j - 1, i]) * (bh[j + 1, i] - bh[j - 1, i])) * (
                        bh_err[j + 1, i] * bh_err[j + 1, i] + bh_err[j - 1, i] * bh_err[j - 1, i]))

    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero.
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the conf_disambig and bitmap arrays.

    i = 0
    for j in range(ny):
        derx_bh[j, i] = ((-3 * bh[j, i]) + (4 * bh[j, i + 1]) - (bh[j, i + 2])) * 0.5

    i = nx - 1
    for j in range(ny):
        derx_bh[j, i] = ((3 * bh[j, i]) + (-4 * bh[j, i - 1]) - (-bh[j, i - 2])) * 0.5

    j = 0
    for i in range(nx):
        dery_bh[j, i] = ((-3 * bh[j, i]) + (4 * bh[j + 1, i]) - (bh[(j + 2), i])) * 0.5

    j = ny - 1
    for i in range(nx):
        dery_bh[j, i] = ((3 * bh[j, i]) + (-4 * bh[j - 1, i]) - (-bh[j - 2, i])) * 0.5

    # Calculate the sum only
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if ((derx_bh[j, i] + dery_bh[j, i]) == 0):
                continue
            if np.isnan(bh[j, i]):
                continue
            if np.isnan(bh[j + 1, i]):
                continue
            if np.isnan(bh[j - 1, i]):
                continue
            if np.isnan(bh[j, i - 1]):
                continue
            if np.isnan(bh[j, i + 1]):
                continue
            if np.isnan(bh_err[j, i]):
                continue
            if np.isnan(derx_bh[j, i]):
                continue
            if np.isnan(dery_bh[j, i]):
                continue
            sum += np.sqrt(derx_bh[j, i] * derx_bh[j, i] + dery_bh[j, i] * dery_bh[j, i])
            err += err_term2[j, i] / (16.0 * (derx_bh[j, i] * derx_bh[j, i] + dery_bh[j, i] * dery_bh[j, i])) + \
                   err_term1[j, i] / (16.0 * (derx_bh[j, i] * derx_bh[j, i] + dery_bh[j, i] * dery_bh[j, i]))
            count_mask += 1

    mean_derivative_bh = (sum) / (count_mask)
    mean_derivative_bh_err = (np.sqrt(err)) / (count_mask)

    return [mean_derivative_bh, mean_derivative_bh_err]


# ===========================================

def computeBzderivative(bz, bz_err, nx, ny, conf_disambig, bitmap):
    """function: computeBzderivative

    This function computes the derivative of the vertical field, or sqrt[(dB_z/dx)^2 + (dB_z/dy)^2].
    The native units in the series hmi.sharp_720s and hmi.sharp_cea_720s are in Gauss/pixel.

    Here are the steps to convert from Gauss/pixel to Gauss/Mm:
    The units of the magnetic field, or dB_z, are in Gauss.
    The units of length, i.e. dx or dy, are in pixels.
    The units of dB_z/dx or dB_z/dy = (Gauss/pix)(pix/arcsec)(arsec/meter)(meter/Mm), or
                                    = (Gauss/pix)(1/cdelt1_arcsec)(RSUN_OBS/RSUN_REF)(1000000)
                                    = Gauss/Mm
    In other words, multiply MEANGBZ by a factor of (1/cdelt1_arcsec)*(RSUN_OBS/RSUN_REF)*(1000000).

    Note that cdelt1_arcsec is defined in the get_data function above.

    """

    count_mask = 0
    sum = 0.0
    err = 0.0

    derx_bz = np.zeros([ny, nx])
    dery_bz = np.zeros([ny, nx])
    err_term1 = np.zeros([ny, nx])
    err_term2 = np.zeros([ny, nx])

    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1, nx - 1):
        for j in range(0, ny):
            derx_bz[j, i] = (bz[j, i + 1] - bz[j, i - 1]) * 0.5
            err_term1[j, i] = (((bz[j, i + 1] - bz[j, i - 1]) * (bz[j, i + 1] - bz[j, i - 1])) * (
                        bz_err[j, i + 1] * bz_err[j, i + 1] + bz_err[j, i - 1] * bz_err[j, i - 1]))

    # brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0, nx):
        for j in range(1, ny - 1):
            dery_bz[j, i] = (bz[j + 1, i] - bz[j - 1, i]) * 0.5
            err_term2[j, i] = (((bz[j + 1, i] - bz[j - 1, i]) * (bz[j + 1, i] - bz[j - 1, i])) * (
                        bz_err[j + 1, i] * bz_err[j + 1, i] + bz_err[j - 1, i] * bz_err[j - 1, i]))

    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero.
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the conf_disambig and bitmap arrays.

    i = 0
    for j in range(ny):
        derx_bz[j, i] = ((-3 * bz[j, i]) + (4 * bz[j, i + 1]) - (bz[j, i + 2])) * 0.5

    i = nx - 1
    for j in range(ny):
        derx_bz[j, i] = ((3 * bz[j, i]) + (-4 * bz[j, i - 1]) - (-bz[j, i - 2])) * 0.5

    j = 0
    for i in range(nx):
        dery_bz[j, i] = ((-3 * bz[j, i]) + (4 * bz[j + 1, i]) - (bz[(j + 2), i])) * 0.5

    j = ny - 1
    for i in range(nx):
        dery_bz[j, i] = ((3 * bz[j, i]) + (-4 * bz[j - 1, i]) - (-bz[j - 2, i])) * 0.5

    # Calculate the sum only
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if ((derx_bz[j, i] + dery_bz[j, i]) == 0):
                continue
            if np.isnan(bz[j, i]):
                continue
            if np.isnan(bz[j + 1, i]):
                continue
            if np.isnan(bz[j - 1, i]):
                continue
            if np.isnan(bz[j, i - 1]):
                continue
            if np.isnan(bz[j, i + 1]):
                continue
            if np.isnan(bz_err[j, i]):
                continue
            if np.isnan(derx_bz[j, i]):
                continue
            if np.isnan(dery_bz[j, i]):
                continue
            sum += np.sqrt(derx_bz[j, i] * derx_bz[j, i] + dery_bz[j, i] * dery_bz[j, i])
            err += err_term2[j, i] / (16.0 * (derx_bz[j, i] * derx_bz[j, i] + dery_bz[j, i] * dery_bz[j, i])) + \
                   err_term1[j, i] / (16.0 * (derx_bz[j, i] * derx_bz[j, i] + dery_bz[j, i] * dery_bz[j, i]))
            count_mask += 1

    mean_derivative_bz = (sum) / (count_mask)
    mean_derivative_bz_err = (np.sqrt(err)) / (count_mask)

    return [mean_derivative_bz, mean_derivative_bz_err]


# ===========================================

def computeJz(bx, by, bx_err, by_err, conf_disambig, bitmap, nx, ny):
    """function: computeJz

    This function computes the z-component of the current.

    In discretized space like data pixels, the current (or curl of B) is calculated as the integration
    of the field Bx and By along the circumference of the data pixel divided by the area of the pixel.

    One form of differencing the curl is expressed as:
    (dx * (Bx(i,j-1)+Bx(i,j)) / 2
    +dy * (By(i+1,j)+By(i,j)) / 2
    -dx * (Bx(i,j+1)+Bx(i,j)) / 2
    -dy * (By(i-1,j)+By(i,j)) / 2) / (dx * dy)

    To change units from Gauss/pixel to mA/m^2 (the units for Jz in Leka and Barnes, 2003),
    one must perform the following unit conversions:
    (Gauss)(1/arcsec)(arcsec/meter)(Newton/Gauss*Ampere*meter)(Ampere^2/Newton)(milliAmpere/Ampere), or
    (Gauss)(1/CDELT1)(RSUN_OBS/RSUN_REF)(1 T / 10^4 Gauss)(1 / 4*PI*10^-7)( 10^3 milliAmpere/Ampere), or
    (Gauss)(1/CDELT1)(RSUN_OBS/RSUN_REF)(0.00010)(1/MUNAUGHT)(1000.),
    where a Tesla is represented as a Newton/Ampere*meter.

    The units of total unsigned vertical current (us_i) are simply in A. In this case, we would have the following:
    (Gauss/pix)(1/CDELT1)(RSUN_OBS/RSUN_REF)(0.00010)(1/MUNAUGHT)(CDELT1)(CDELT1)(RSUN_REF/RSUN_OBS)(RSUN_REF/RSUN_OBS)
    = (Gauss/pix)(0.00010)(1/MUNAUGHT)(CDELT1)(RSUN_REF/RSUN_OBS)
    """

    derx = np.zeros([ny, nx])
    dery = np.zeros([ny, nx])
    err_term1 = np.zeros([ny, nx])
    err_term2 = np.zeros([ny, nx])
    jz = np.zeros([ny, nx])
    jz_err = np.zeros([ny, nx])

    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1, nx - 1):
        for j in range(0, ny):
            derx[j, i] = (by[j, i + 1] - by[j, i - 1]) * 0.5
            err_term1[j, i] = ((by_err[j, i + 1] * by_err[j, i + 1]) + (by_err[j, i - 1] * by_err[j, i - 1]))

    # brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0, nx):
        for j in range(1, ny - 1):
            dery[j, i] = (bx[j + 1, i] - bx[j - 1, i]) * 0.5
            err_term2[j, i] = ((bx_err[j + 1, i] * bx_err[j + 1, i]) + (bx_err[j - 1, i] * bx_err[j - 1, i]))

    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero.
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the conf_disambig and bitmap arrays.

    i = 0
    for j in range(ny):
        derx[j, i] = ((-3 * by[j, i]) + (4 * by[j, i + 1]) - (by[j, i + 2])) * 0.5

    i = nx - 1
    for j in range(ny):
        derx[j, i] = ((3 * by[j, i]) + (-4 * by[j, i - 1]) - (-by[j, i - 2])) * 0.5

    j = 0
    for i in range(nx):
        dery[j, i] = ((-3 * bx[j, i]) + (4 * bx[j + 1, i]) - (bx[(j + 2), i])) * 0.5

    j = ny - 1
    for i in range(nx):
        dery[j, i] = ((3 * bx[j, i]) + (-4 * bx[j - 1, i]) - (-bx[j - 2, i])) * 0.5

    # Calculate the sum only
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            jz[j, i] = (derx[j, i] - dery[j, i])
            jz_err[j, i] = 0.5 * np.sqrt(err_term1[j, i] + err_term2[j, i])

    return [jz, jz_err, derx, dery]


# ===========================================

def computeJzmoments(jz, jz_err, derx, dery, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec,
                     munaught):
    """function: computeJzmoments

    This function computes moments of the vertical current.
    The mean vertical current density is in units of mA/m^2.
    The total unsigned vertical current is in units of Amperes.
    """

    count_mask = 0
    curl = 0.0
    err = 0.0
    us_i = 0.0

    # Calculate the sum only
    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(jz[j, i]):
                continue
            if np.isnan(derx[j, i]):
                continue
            if np.isnan(dery[j, i]):
                continue
            curl += (jz[j, i]) * (1 / cdelt1_arcsec) * (rsun_obs / rsun_ref) * (0.00010) * (1 / munaught) * (1000.)
            us_i += abs(jz[j, i]) * (cdelt1_arcsec / 1) * (rsun_ref / rsun_obs) * (0.00010) * (1 / munaught)
            err += (jz_err[j, i] * jz_err[j, i])
            count_mask += 1

    mean_jz = curl / (count_mask)
    mean_jz_err = (np.sqrt(err) / count_mask) * (
                (1 / cdelt1_arcsec) * (rsun_obs / rsun_ref) * (0.00010) * (1 / munaught) * (1000.))

    us_i = (us_i)
    us_i_err = (np.sqrt(err)) * ((cdelt1_arcsec / 1) * (rsun_ref / rsun_obs) * (0.00010) * (1 / munaught))

    return [mean_jz, mean_jz_err, us_i, us_i_err]


# ===========================================

def computeAlpha(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec):
    """function: computeAlpha

    This function computes the twist parameter.

    The twist parameter, alpha, is defined as alpha = Jz/Bz. In this case, the calculation for alpha is weighted by Bz:

    numerator   = sum of all Jz*Bz
    denominator = sum of Bz*Bz
    alpha       = numerator/denominator

    The units of alpha are in 1/Mm
    The units of Jz are in Gauss/pix; the units of Bz are in Gauss.

    Therefore, the units of Jz/Bz = (Gauss/pix)(1/Gauss)(pix/arcsec)(arsec/meter)(meter/Mm), or
    = (Gauss/pix)(1/Gauss)(1/CDELT1)(RSUN_OBS/RSUN_REF)(10^6)
    = 1/Mm
    """

    alpha_total = 0.0
    C = ((1 / cdelt1_arcsec) * (rsun_obs / rsun_ref) * (1000000.))
    total = 0.0
    A = 0.0
    B = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(jz[j, i]):
                continue
            if np.isnan(bz[j, i]):
                continue
            if (jz[j, i] == 0):
                continue
            if (bz[j, i] == 0):
                continue
            A += jz[j, i] * bz[j, i]
            B += bz[j, i] * bz[j, i]

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(jz[j, i]):
                continue
            if np.isnan(bz[j, i]):
                continue
            if (jz[j, i] == 0):
                continue
            if (bz[j, i] == 0):
                continue
            total += bz[j, i] * bz[j, i] * jz_err[j, i] * jz_err[j, i] + (jz[j, i] - 2 * bz[j, i] * A / B) * (
                        jz[j, i] - 2 * bz[j, i] * A / B) * bz_err[j, i] * bz_err[j, i]

    # Determine the absolute value of alpha. The units for alpha are 1/Mm
    alpha_total = ((A / B) * C)
    mean_alpha = alpha_total
    mean_alpha_err = (C / B) * (np.sqrt(total))

    return [mean_alpha, mean_alpha_err]


# ===========================================

def computeHelicity(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec):
    """function: computeHelicity

    This function computes a proxy for the current helicity and various moments.

    The current helicity is defined as Bz*Jz and the units are G^2 / m
    The units of Jz are in G/pix; the units of Bz are in G.
    Therefore, the units of Bz*Jz = (Gauss)*(Gauss/pix) = (Gauss^2/pix)(pix/arcsec)(arcsec/meter)
    = (Gauss^2/pix)(1/CDELT1)(RSUN_OBS/RSUN_REF)
    =  G^2 / m.
    """

    count_mask = 0.0
    sum = 0.0
    sum2 = 0.0
    err = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(jz[j, i]):
                continue
            if np.isnan(bz[j, i]):
                continue
            if (jz[j, i] == 0):
                continue
            if (bz[j, i] == 0):
                continue
            if np.isnan(jz_err[j, i]):
                continue
            if np.isnan(bz_err[j, i]):
                continue
            sum += (jz[j, i] * bz[j, i]) * (1 / cdelt1_arcsec) * (
                        rsun_obs / rsun_ref)  # contributes to MEANJZH and ABSNJZH
            sum2 += abs(jz[j, i] * bz[j, i]) * (1 / cdelt1_arcsec) * (rsun_obs / rsun_ref)  # contributes to TOTUSJH
            err += (jz_err[j, i] * jz_err[j, i] * bz[j, i] * bz[j, i]) + (
                        bz_err[j, i] * bz_err[j, i] * jz[j, i] * jz[j, i])
            count_mask += 1

    mean_ih = sum / count_mask  # Units are G^2 / m ; keyword is MEANJZH
    total_us_ih = sum2  # Units are G^2 / m ; keyword is TOTUSJH
    total_abs_ih = abs(sum)  # Units are G^2 / m ; keyword is ABSNJZH
    mean_ih_err = (np.sqrt(err) / count_mask) * (1 / cdelt1_arcsec) * (
                rsun_obs / rsun_ref)  # error in the quantity MEANJZH
    total_us_ih_err = (np.sqrt(err)) * (1 / cdelt1_arcsec) * (rsun_obs / rsun_ref)  # error in the quantity TOTUSJH
    total_abs_ih_err = (np.sqrt(err)) * (1 / cdelt1_arcsec) * (rsun_obs / rsun_ref)  # error in the quantity ABSNJZH

    return [mean_ih, mean_ih_err, total_us_ih, total_us_ih_err, total_abs_ih, total_abs_ih_err]


# ===========================================

def computeSumAbsPerPolarity(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec,
                             munaught):
    """function: computeSumAbsPerPolarity

    This function computes the sum of the absolute value of the current per polarity. It is defined as follows:

    The sum of the absolute value per polarity is defined as the following:
    abs(sum(jz gt 0)) + abs(sum(jz lt 0)) and the units are in Amperes per arcsecond.
    The units of jz are in G/pix. In this case, we would have the following:
    Jz = (Gauss/pix)(1/CDELT1)(0.00010)(1/MUNAUGHT)(RSUN_REF/RSUN_OBS)(RSUN_REF/RSUN_OBS)(RSUN_OBS/RSUN_REF),
       = (Gauss/pix)(1/CDELT1)(0.00010)(1/MUNAUGHT)(RSUN_REF/RSUN_OBS)

    The error in this quantity is the same as the error in the mean vertical current.
    """

    count_mask = 0.0
    sum1 = 0.0
    sum2 = 0.0
    err = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(jz[j, i]):
                continue
            if np.isnan(bz[j, i]):
                continue
            if (bz[j, i] > 0):
                sum1 += (jz[j, i]) * (1 / cdelt1_arcsec) * (0.00010) * (1 / munaught) * (rsun_ref / rsun_obs)
            if (bz[j, i] <= 0):
                sum2 += (jz[j, i]) * (1 / cdelt1_arcsec) * (0.00010) * (1 / munaught) * (rsun_ref / rsun_obs)
            err += (jz_err[j, i] * jz_err[j, i])
            count_mask += 1

    totaljz = abs(sum1) + abs(sum2)
    totaljz_err = np.sqrt(err) * (1 / cdelt1_arcsec) * abs((0.00010) * (1 / munaught) * (rsun_ref / rsun_obs))

    return [totaljz, totaljz_err]


# ===========================================

def computeFreeEnergy(bx_err, by_err, bx, by, bpx, bpy, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec, conf_disambig,
                      bitmap):
    """
    function: computeFreeEnergy

    This function computes the mean photospheric excess magnetic energy and total photospheric excess magnetic energy density.

    The units for magnetic energy density in cgs are ergs per cubic centimeter. The formula B^2/8*PI integrated over all space, dV
    automatically yields erg per cubic centimeter for an input B in Gauss. Note that the 8*PI can come out of the integral; thus,
    the integral is over B^2 dV and the 8*PI is divided at the end.

    Total magnetic energy is the magnetic energy density times dA, or the area, and the units are thus ergs/cm. To convert
    ergs per centimeter cubed to ergs per centimeter, simply multiply by the area per pixel in cm:
    erg/cm^3*(CDELT1^2)*(RSUN_REF/RSUN_OBS ^2)*(100.^2)
    = erg/cm(1/pix^2)
    """
    count_mask = 0.0
    sum = 0.0
    sum1 = 0.0
    err = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(bx[j, i]):
                continue
            if np.isnan(by[j, i]):
                continue
            sum += (((bx[j, i] - bpx[j, i]) * (bx[j, i] - bpx[j, i])) + (
                        (by[j, i] - bpy[j, i]) * (by[j, i] - bpy[j, i]))) * (
                               cdelt1_arcsec * cdelt1_arcsec * (rsun_ref / rsun_obs) * (
                                   rsun_ref / rsun_obs) * 100.0 * 100.0)
            sum1 += (((bx[j, i] - bpx[j, i]) * (bx[j, i] - bpx[j, i])) + (
                        (by[j, i] - bpy[j, i]) * (by[j, i] - bpy[j, i])))
            err += 4.0 * (bx[j, i] - bpx[j, i]) * (bx[j, i] - bpx[j, i]) * (bx_err[j, i] * bx_err[j, i]) + 4.0 * (
                        by[j, i] - bpy[j, i]) * (by[j, i] - bpy[j, i]) * (by_err[j, i] * by_err[j, i])
            count_mask += 1

    # Units of meanpotptr are ergs per centimeter
    meanpot = (sum1) / (count_mask * 8. * np.pi)
    meanpot_err = (np.sqrt(err)) / (count_mask * 8. * np.pi)

    # Units of sum are ergs/cm^3, units of factor are cm^2/pix^2; therefore, units of totpotptr are ergs per centimeter
    totpot = (sum) / (8. * np.pi)
    totpot_err = (np.sqrt(err)) * abs(
        cdelt1_arcsec * cdelt1_arcsec * (rsun_ref / rsun_obs) * (rsun_ref / rsun_obs) * 100.0 * 100.0 * (
                    1 / (8. * np.pi)))

    return [meanpot, meanpot_err, totpot, totpot_err]


# ===========================================

def computeShearAngle(bx_err, by_err, bz_err, bx, by, bz, bpx, bpy, nx, ny, conf_disambig, bitmap):
    """
    function: computeShearAngle

    This function computes the shear angle, or the angle between the potential field vector and the observed field vector, in degrees.
    """

    count_mask = 0.0
    count = 0.0
    dotproduct = 0.0
    magnitude_potential = 0.0
    magnitude_vector = 0.0
    sumsum = 0.0
    shear_angle = 0.0
    denominator = 0.0
    term1 = 0.0
    sumsum = 0.0
    err = 0.0
    part1 = 0.0
    part2 = 0.0
    part3 = 0.0

    for j in range(ny):
        for i in range(nx):
            if (conf_disambig[j, i] < 70 or bitmap[j, i] < 30):
                continue
            if np.isnan(bx[j, i]):
                continue
            if np.isnan(by[j, i]):
                continue
            if np.isnan(bz[j, i]):
                continue
            if np.isnan(bpx[j, i]):
                continue
            if np.isnan(bpy[j, i]):
                continue
            if np.isnan(bx_err[j, i]):
                continue
            if np.isnan(by_err[j, i]):
                continue
            if np.isnan(bz_err[j, i]):
                continue
            # for the values
            dotproduct = (bpx[j, i]) * (bx[j, i]) + (bpy[j, i]) * (by[j, i]) + (bz[j, i]) * (bz[j, i])
            magnitude_potential = np.sqrt((bpx[j, i] * bpx[j, i]) + (bpy[j, i] * bpy[j, i]) + (bz[j, i] * bz[j, i]))
            magnitude_vector = np.sqrt((bx[j, i] * bx[j, i]) + (by[j, i] * by[j, i]) + (bz[j, i] * bz[j, i]))
            shear_angle = math.acos(dotproduct / (magnitude_potential * magnitude_vector)) * (180. / np.pi)
            sumsum += shear_angle
            count += 1
            # for the error analysis
            term1 = bx[j, i] * by[j, i] * bpy[j, i] - by[j, i] * by[j, i] * bpx[j, i] + bz[j, i] * bx[j, i] * bz[j, i] - \
                    bz[j, i] * bz[j, i] * bpx[j, i]
            # term2 = bx[j,i]*bx[j,i]*bpy[j,i] - bx[j,i]*by[j,i]*bpx[j,i] + bx[j,i]*bz[j,i]*bpy[j,i] - bz[j,i]*by[j,i]*bz[j,i]
            # term3 = bx[j,i]*bx[j,i]*bz[j,i] - bx[j,i]*bz[j,i]*bpx[j,i] + by[j,i]*by[j,i]*bz[j,i] - by[j,i]*bz[j,i]*bpy[j,i]
            part1 = bx[j, i] * bx[j, i] + by[j, i] * by[j, i] + bz[j, i] * bz[j, i]
            part2 = bpx[j, i] * bpx[j, i] + bpy[j, i] * bpy[j, i] + bz[j, i] * bz[j, i]
            part3 = bx[j, i] * bpx[j, i] + by[j, i] * bpy[j, i] + bz[j, i] * bz[j, i]
            denominator = part1 * part1 * part1 * part2 * (1.0 - ((part3 * part3) / (part1 * part2)))
            err = (term1 * term1 * bx_err[j, i] * bx_err[j, i]) / (denominator) + (
                        term1 * term1 * bx_err[j, i] * bx_err[j, i]) / (denominator) + (
                              term1 * term1 * bx_err[j, i] * bx_err[j, i]) / (denominator)
            if (shear_angle > 45):
                count_mask += 1

    # For mean 3D shear angle, area with shear greater than 45
    meanshear_angle = (sumsum) / (count)
    meanshear_angle_err = (np.sqrt(err) / count_mask) * (180. / np.pi)

    # The area here is a fractional area -- the % of the total area. This has no error associated with it.
    area_w_shear_gt_45 = (count_mask / (count)) * (100.0)

    return [meanshear_angle, meanshear_angle_err, area_w_shear_gt_45]


# ===========================================

def computeR(los, los_err, nx, ny, cdelt1_arcsec):
    """
    function: computeR

    This function computes R, or the log of the gradient-weighted neutral line length.
    So the output is unitless.

    This function also computes the error in R. The general formula for an error of a
    function is ERR(f(x)) = d/dx f(x) * ERR(x). Thus
                ERR(R) = d/dx (log_10 R) * ERR(x)
                       = [1/ln(10)] * [1/x] * ERR(x)
                       = ERR(x) / ln(10)*x
    """

    sum = 0.0
    err = 0.0
    sigma = 10.0 / 2.3548
    scale = int(round(2.0 / cdelt1_arcsec))

    # =============== [STEP 1] ===============
    # bin the line-of-sight magnetogram down by a factor of scale
    rim = block_reduce(los, block_size=(scale, scale), func=np.mean)

    # =============== [STEP 2] ===============
    # identify positive and negative pixels greater than +/- 150 gauss
    # and label those pixels with a 1.0 in arrays p1p0 and p1n0

    nx1 = rim.shape[1]
    ny1 = rim.shape[0]
    p1p0 = np.zeros([ny1, nx1])
    p1n0 = np.zeros([ny1, nx1])

    for j in range(ny1):
        for i in range(nx1):
            if (rim[j, i] > 150):
                p1p0[j, i] = 1.0
            else:
                p1p0[j, i] = 0.0
            if (rim[j, i] < -150):
                p1n0[j, i] = 1.0
            else:
                p1n0[j, i] = 0.0

    # =============== [STEP 3] ===============
    # smooth each of the negative and positive pixel bitmaps by convolving with a boxcar

    # set up the convolution kernel
    boxcar_kernel = np.zeros([ny1, nx1])
    midpoint_ny1 = int(round(ny1 / 2))
    midpoint_nx1 = int(round(nx1 / 2))

    for j in range(midpoint_ny1, midpoint_ny1 + 3):
        for i in range(midpoint_nx1, midpoint_nx1 + 3):
            boxcar_kernel[j, i] = 0.1111

    p1p = scipy.ndimage.convolve(p1p0, boxcar_kernel)
    p1n = scipy.ndimage.convolve(p1n0, boxcar_kernel)

    # =============== [STEP 4] ===============
    # find the pixels for which p1p and p1n are both equal to 1.
    # this defines the polarity inversion line

    p1 = np.zeros([ny1, nx1])
    for j in range(ny1):
        for i in range(nx1):
            if ((p1p[j, i] > 0.0) and (p1n[j, i] > 0.0)):
                p1[j, i] = 1.0
            else:
                p1[j, i] = 0.0

    # =============== [STEP 5] ===============
    # convolve the polarity inversion line map with a gaussian
    # to identify the region near the plarity inversion line
    # the resultant array is called pmap

    pmap = scipy.ndimage.gaussian_filter(p1, sigma, order=0)

    # =============== [STEP 6] ===============
    # the R parameter is calculated

    for j in range(ny1):
        for i in range(nx1):
            if np.isnan(pmap[j, i]):
                continue
            if np.isnan(rim[j, i]):
                continue
            sum += pmap[j, i] * abs(rim[j, i])
            err += pmap[j, i] * abs(los_err[j, i])

    if (sum < 1.0):
        Rparam = 0.0
        Rparam_err = 0.0
    else:
        Rparam = math.log10(sum)
        Rparam_err = err / (math.log(10) * sum)  # note that math.log is a natural log by default

    return [Rparam, Rparam_err]


# ===========================================

def computeLOSderivative(los, los_err, nx, ny, bitmap, rsun_ref, rsun_obs, cdelt1_arcsec):
    """function: computeLOSderivative

    This function computes the derivative of the line-of-sight field, or sqrt[(dB_los/dx)^2 + (dB_los/dy)^2].
    The native units in the series hmi.sharp_720s and hmi.sharp_cea_720s are in Gauss/pixel.

    Here are the steps to convert from Gauss/pixel to Gauss/Mm:
    The units of the magnetic field, or dB_los, are in Gauss.
    The units of length, i.e. dx or dy, are in pixels.
    The units of dB_los/dx or dB_los/dy = (Gauss/pix)(pix/arcsec)(arsec/meter)(meter/Mm), or
                                        = (Gauss/pix)(1/cdelt1_arcsec)(RSUN_OBS/RSUN_REF)(1000000)
                                        = Gauss/Mm
    In other words, multiply MEANGBL by a factor of (1/cdelt1_arcsec)*(RSUN_OBS/RSUN_REF)*(1000000).

    Note that cdelt1_arcsec is defined in the get_data function above.
    """

    count_mask = 0
    sum = 0.0
    err = 0.0
    unitconstant = (1 / cdelt1_arcsec) * (rsun_obs / rsun_ref) * (10e6)

    derx_blos = np.zeros([ny, nx])
    dery_blos = np.zeros([ny, nx])
    err_term1 = np.zeros([ny, nx])
    err_term2 = np.zeros([ny, nx])

    # brute force method of calculating the derivative d/dx (no consideration for edges)
    for i in range(1, nx - 1):
        for j in range(0, ny):
            derx_blos[j, i] = (los[j, i + 1] - los[j, i - 1]) * 0.5
            err_term1[j, i] = (((los[j, i + 1] - los[j, i - 1]) * (los[j, i + 1] - los[j, i - 1])) * (
                        los_err[j, i + 1] * los_err[j, i + 1] + los_err[j, i - 1] * los_err[j, i - 1]))

    # brute force method of calculating the derivative d/dy (no consideration for edges) */
    for i in range(0, nx):
        for j in range(1, ny - 1):
            dery_blos[j, i] = (los[j + 1, i] - los[j - 1, i]) * 0.5
            err_term2[j, i] = (((los[j + 1, i] - los[j - 1, i]) * (los[j + 1, i] - los[j - 1, i])) * (
                        los_err[j + 1, i] * los_err[j + 1, i] + los_err[j - 1, i] * los_err[j - 1, i]))

    # consider the edges for the arrays that contribute to the variable "sum" in the computation below.
    # ignore the edges for the error terms as those arrays have been initialized to zero.
    # this is okay because the error term will ultimately not include the edge pixels as they are selected out by the conf_disambig and bitmap arrays.

    i = 0
    for j in range(ny):
        derx_blos[j, i] = ((-3 * los[j, i]) + (4 * los[j, i + 1]) - (los[j, i + 2])) * 0.5

    i = nx - 1
    for j in range(ny):
        derx_blos[j, i] = ((3 * los[j, i]) + (-4 * los[j, i - 1]) - (-los[j, i - 2])) * 0.5

    j = 0
    for i in range(nx):
        dery_blos[j, i] = ((-3 * los[j, i]) + (4 * los[j + 1, i]) - (los[(j + 2), i])) * 0.5

    j = ny - 1
    for i in range(nx):
        dery_blos[j, i] = ((3 * los[j, i]) + (-4 * los[j - 1, i]) - (-los[j - 2, i])) * 0.5

    # Calculate the sum only
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            if (bitmap[j, i] < 30):
                continue
            if np.isnan(los[j, i]):
                continue
            if np.isnan(los[j + 1, i]):
                continue
            if np.isnan(los[j - 1, i]):
                continue
            if np.isnan(los[j, i - 1]):
                continue
            if np.isnan(los[j, i + 1]):
                continue
            if np.isnan(derx_blos[j, i]):
                continue
            if np.isnan(dery_blos[j, i]):
                continue
            sum += np.sqrt(derx_blos[j, i] * derx_blos[j, i] + dery_blos[j, i] * dery_blos[j, i])
            denominator_1 = 16.0 * (derx_blos[j, i] * derx_blos[j, i] + dery_blos[j, i] * dery_blos[j, i])
            denominator_2 = 16.0 * (derx_blos[j, i] * derx_blos[j, i] + dery_blos[j, i] * dery_blos[j, i])
            if np.isnan(denominator_1):
                continue
            if np.isnan(denominator_2):
                continue
            if denominator_1 == 0:
                continue
            if denominator_2 == 0:
                continue
            err += (err_term2[j, i] / denominator_1) + (err_term1[j, i] / denominator_2)
            count_mask += 1

    if count_mask == 0:
        mean_derivative_blos = 0.0
        mean_derivative_blos_err = 0.0
    else:
        mean_derivative_blos = (sum) / (count_mask)
        mean_derivative_blos_err = (np.sqrt(err)) / (count_mask)

    return [mean_derivative_blos, mean_derivative_blos_err]


# ===========================================

def compute_abs_flux_los(los, los_err, bitmap, nx, ny, rsun_ref, rsun_obs, cdelt1_arcsec):
    """function: compute_abs_flux_los

    This function computes the total unsigned flux, on the line-of-sight field, in units of G/cm^2.
    It also returns the number of pixels used in this calculation in the keyword CMASK.

    To compute the unsigned flux, we simply calculate
       flux = surface integral [(vector Blos) dot (normal vector)],
            = surface integral [(magnitude Blos)*(magnitude normal)*(cos theta)].

    However, since the field is radial, we will assume cos theta = 1.
    Therefore, the pixels only need to be corrected for the projection.

    To convert G to G*cm^2, simply multiply by the number of square centimeters per pixel:
       (Gauss/pix^2)(CDELT1)^2(RSUN_REF/RSUN_OBS)^2(100.cm/m)^2
       =Gauss*cm^2
    """

    count_mask = 0
    sum = 0.0
    err = 0.0

    for j in range(ny):
        for i in range(nx):
            if (bitmap[j, i] < 30):
                continue
            if np.isnan(los[j, i]):
                continue
            sum += abs(los[j, i])
            err += los_err[j, i] * los_err[j, i]
            count_mask += 1

    mean_vf = sum * cdelt1_arcsec * cdelt1_arcsec * (rsun_ref / rsun_obs) * (rsun_ref / rsun_obs) * 100.0 * 100.0
    mean_vf_err = (np.sqrt(err)) * abs(
        cdelt1_arcsec * cdelt1_arcsec * (rsun_ref / rsun_obs) * (rsun_ref / rsun_obs) * 100.0 * 100.0)

    return [mean_vf, mean_vf_err, count_mask]


# ===========================================

def greenpot(bz, nx, ny):
    """
    function: greenpot

    This function extrapolates the potential magnetic field using Green's functions.
    The underlying assuption of a potential field is that it is Maxwell-stress free.
    The monopole depth is 0.01 pixels.
    """
    print('Calculating the potential field. This takes a minute.')

    nnx = nx
    nny = ny

    # define the monopole depth, dz
    dz = 0.001

    # malloc some arrays
    pfpot = np.zeros([nny, nnx])
    rdist = np.zeros([nny, nnx])
    bztmp = np.zeros([nny, nnx])
    bxp = np.zeros([nny, nnx])
    byp = np.zeros([nny, nnx])

    # substitute zeros for nans in bz
    for iny in range(nny):
        for inx in range(nnx):
            if np.isnan(bz[iny, inx]):
                bztmp[iny, inx] = 0.0
            else:
                bztmp[iny, inx] = bz[iny, inx]

    rdd = 0.0
    rdd1 = 0.0
    rdd2 = 0.0
    for iny in range(nny):
        for inx in range(nnx):
            rdd1 = float(inx)
            rdd2 = float(iny)
            rdd = rdd1 * rdd1 + rdd2 * rdd2 + dz * dz
            rdist[iny, inx] = 1.0 / (np.sqrt(rdd))

    iwindow = 0
    if (nnx > nny):
        iwindow = nnx
    else:
        iwindow = nny

    rwindow = float(iwindow)
    rwindow = rwindow * rwindow + 0.01  # must be square
    rwindow = 1.0e2  # limit the window size to be 10.
    rwindow = np.sqrt(rwindow)
    iwindow = int(rwindow)

    for iny in range(nny):
        for inx in range(nnx):
            if np.isnan(bz[iny, inx]):
                pfpot[iny, inx] = 0.0
            else:
                sum = 0.0
                j2s = iny - iwindow
                j2e = iny + iwindow
                if (j2s < 0):
                    j2s = 0
                if (j2e > nny):
                    j2e = nny
                i2s = inx - iwindow
                i2e = inx + iwindow
                if (i2s < 0):
                    i2s = 0
                if (i2e > nnx):
                    i2e = nnx
                for j2 in range(j2s, j2e):
                    for i2 in range(i2s, i2e):
                        val1 = bztmp[j2, i2]
                        di = abs(i2 - inx)
                        dj = abs(j2 - iny)
                        sum = sum + val1 * rdist[dj, di] * dz
                pfpot[iny, inx] = sum

    for iny in range(1, nny - 1):
        for inx in range(1, nnx - 1):
            bxp[iny, inx] = -(pfpot[iny, inx + 1] - pfpot[iny, inx - 1]) * 0.5
            byp[iny, inx] = -(pfpot[iny + 1, inx] - pfpot[iny - 1, inx]) * 0.5

    return [bxp, byp]

# ===========================================