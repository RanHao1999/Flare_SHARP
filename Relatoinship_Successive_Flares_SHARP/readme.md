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

3. Go to 'SHARP_original(pilmaksed, original_calculated).py', which are written for compute SHARP parameters with original B data  or PIL-masked B data.
   Original: SHARPs are directly read from the .fits file, original_calculated: read the B data and calculate the SHARPs in the whole AR, pilmaksed: read the B data and calculate the SHARPs in the PIL area.
   
   At this step, we output .csv files each carries the SHARP parameters at one time.

4. Run 'csv_read.py', we output two .csv files, which contain all SHARP parameters for one sort of selected data.

5. Run 'Flares', one obtains the final result.

The results of the programme is published in ApJ 937:43.

Discussion and suggestions about the codes are warmly welcomed!

#####################################################

Enjoy!
