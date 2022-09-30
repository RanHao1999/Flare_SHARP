#####################################################

This is the Readme file in this program.

__author__ = @RanHao

The functions for calculating SHARP parameters are adopted from Prof. Monica Bobra (@mbobra).

Steps for using the codes:

0. Run 'create_folders.py'. The code creates the proper folders.

1. Download HMI magnetic field datas into the 'data/JSOC' folder, which should include 'Br.fits', 'Bp.fits', 'Bt.fits',
   'Br_err.fits', 'Bp_err.fits', 'Bt_err.fits', 'bitmap.fits', 'conf_disambig.fits', 'magnetogram.fits'

2. Run 'Locat_PIL.py', which is written for locating the polarity inversion line area.
   Two parameters in the code could be modified: 'sigma' in 'gauss_filter' in line 191 and line 192, max_workers in line 238.
   The first one adjusts the width of the gaussian kernel. The second one determines the number of cores one wants to use in his computer.
   At this step, we output 'Br.jpg', 'bitmap_positive.jpg', 'bitmap_negative.jpg', 'coord_pil.txt' and 'map_pil.jpg' in the corresponding folders in 'res/...'
   If needed, one can run 'Make_video.py', by which we output the .avi videos of the evolution of 'Br', 'map_pil', 'bimap_positive' and 'bitmap_negative'.

3. Go to 'SHARP_with_PIL.py', which is written for compute SHARP parameters with both original B data and PIL-masked B data.
   Only two parameters needs to be modified: 'mask or not' in line 1262 and 'max_workers' in line 1256.
   if 'mask or not' == 'no', then the results are from original B data (the whole AR), else from PIL-masked data.
   At this step, we output .csv files each carries the SHARP parameters at one time.

4. Run 'csv_read.py', we output two .csv files, which contain all SHARP parameters for one sort of selected data.

5. Run 'Flares', one obtains the final result.

The content of the programme is published in ApJ 937:43.
Discussion and suggestions about the codes are warmly welcomed!
#####################################################

Enjoy!
