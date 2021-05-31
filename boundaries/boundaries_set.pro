pro boundaries_set

m = 5
date_lst = ['20170902','20140328','20120213','20170906','20120327','20110214','20110215']
time_lst = ['151200','171200','193600','202400','173600','233600','033600']
lst_use = [date_lst[m],'_',time_lst[m]]
xrange = [0.0,488.0]
yrange = [-330.0,-130.0]

file_field = 'F:\3D-magnetic-con\3d_myrun\01data\'+strjoin(lst_use)+'\hmi.B_720s.'+strjoin(lst_use)+'_TAI.field.fits'
file_inclination = 'F:\3D-magnetic-con\3d_myrun\01data\'+strjoin(lst_use)+'\hmi.B_720s.'+strjoin(lst_use)+'_TAI.inclination.fits'
file_azimuth = 'F:\3D-magnetic-con\3d_myrun\01data\'+strjoin(lst_use)+'\hmi.B_720s.'+strjoin(lst_use)+'_TAI.azimuth.fits'
file_disambig = 'F:\3D-magnetic-con\3d_myrun\01data\'+strjoin(lst_use)+'\hmi.B_720s.'+strjoin(lst_use)+'_TAI.disambig.fits'
read_sdo,file_field,index1,data1,/noshell
index2map,index1,data1,mapf

xc = 0
yc = 0
nx = 0
ny = 0
plot_hmi_func,file_field,file_inclination,file_azimuth,file_disambig,date_lst[m],time_lst[m],xc,yc
xc = xc
yc = yc
projection_func,date_lst[m],time_lst[m]
b_plot_submap_func,date_lst[m],time_lst[m],xrange,yrange
creb_lv3_func,date_lst[m],time_lst[m],nx,ny
nx = nx & ny = ny
xc1 = xc & xc2 = xc & yc1 = yc & yc2 = yc & nx1 = nx
nx2 = nx & ny1 = ny & ny2 = ny
;nx = ceil(xrange[1]-xrange[0])/2.0
;ny = ceil(xrange[1]-xrange[0])/2.0
;read_boundary_potential_func,date_lst[m],time_lst[m],xc1,yc1,nx1,ny1
;print,'step1',xc,yc
read_boundary_nlfff_func,date_lst[m],time_lst[m],xc2,yc2,nx2,ny2
;make_directory,date_lst[m],time_lst[m]

END