Codes specially written for calculating the boundary conditions of NLFFF and potential field of a solar active region.

boundaries_set.pro is the main procedure, codes in the enviroments branch is the basic functions which should be compiled before use, so are the codes in 'functions'.

The name of data file should be "yyyymmdd_hhmmss", such as the example in 01data.

There are three parameteres to be modified:
1. m: choose the date of the data.
2. xrange and yrange in boundaries_set.pro, some restrictions should be considered. Details refer to @njuguoyang.
3. xrange and yrange in plot_hmi_func.pro, determine the subplot area.