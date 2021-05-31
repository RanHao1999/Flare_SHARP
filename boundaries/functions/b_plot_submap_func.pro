pro b_plot_submap_func,date1,time1,xrange,yrange

time_lst = [date1,'_',time1]
filename = 'F:\3D-magnetic-con\3d_myrun\02ambiguity_hmi\05projection\'+strjoin(time_lst)

restore,filename+'\'+'mapbxyz_m.sav',/ver
;xrange = [150.0,638.0]
;yrange = [-300.0,-66.5] ;set properly to fullfill the requirement ,res = nx + 4
res = [0.5,0.5]
sub_map,mapbx,smapbx,xrange=xrange,yrange=yrange
dsmapbx = drot_map(smapbx,time = smapbx.time,resolution=res,/preserve_area)
smapbx = dsmapbx
sub_map,mapby,smapby,xrange=xrange,yrange=yrange
dsmapby = drot_map(smapby,time = smapby.time,resolution=res,/preserve_area)
smapby = dsmapby
sub_map,mapbz,smapbz,xrange=xrange,yrange=yrange
dsmapbz = drot_map(smapbz,time = smapbz.time,resolution=res,/preserve_area)
smapbz = dsmapbz
wdef,1,1200,1200
!p.background=255
smapbz.id='SDO/HMI Heliographic Bz'

plot_map,smapbz,color=0,bcolor=0,charsize=2.0,dmax=800,dmin=-800,xtitle='X (arcsec)',ytitle='Y (arcsec)'
write_png,filename+'\'+'sbz_m.png',tvrd()
smapbz.id='Heliographic Bxyz'
print,'Flux balance coefficient:',total(smapbz.data)/total(abs(smapbz.data))
plot_map,smapbz,color=0,bcolor=0,charsize=2.0,dmax=800,dmin=-800,xtitle='X (arcsec)',ytitle='Y (arcsec)'
plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=150,scale=0.01,iskip=15,jskip=15,$
  v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/Noaxis
!p.background=0
write_png,filename+'\'+'sbxyz_m.png',tvrd()

save,smapbx,smapby,smapbz,filename = filename+'\'+'smapbxyz_m.sav'

END