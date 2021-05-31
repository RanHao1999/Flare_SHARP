pro plot_hmi_func,file_field,file_inclination,file_azimuth,file_disambig,date1,time1,xc,yc

time_lst = [date1,'_',time1]
spawn,'mkdir '+'F:\3D-magnetic-con\3d_myrun\02ambiguity_hmi\01plot_hmi\'+strjoin(time_lst)
filename = 'F:\3D-magnetic-con\3d_myrun\02ambiguity_hmi\01plot_hmi\'+strjoin(time_lst)
spawn,'mkdir '+filename

read_sdo,file_field,index1,data1,/noshell
index2map,index1,data1,mapb
read_sdo,file_inclination,index2,data2,/noshell
index2map,index2,data2,mapi
read_sdo,file_azimuth,index3,data3,/noshell
read_sdo,file_disambig,index4,data4,/noshell

ns = size(data4)
for j=0,ns[2]-1 do begin
  for i=0,ns[1]-1 do begin 
    if (data4[i,j] gt 3) then begin
      data3[i,j]=data3[i,j]+180.0
    endif
  endfor
endfor
index2map,index3,data3,mapa

mapbz = mapb
mapbz.data = rotate(mapb.data*cos(mapi.data*!dtor),2)
mapbz.roll_angle = 0.0
extreme_gap_ud = find_extreme_up(mapbz.data) - find_extreme_down(mapbz.data)
extreme_gap_lr = find_extreme_left(mapbz.data) - find_extreme_right(mapbz.data)
print,extreme_gap_ud,extreme_gap_lr

dx0 = 0.5*(extreme_gap_ud)
dy0 = 0.5*extreme_gap_lr
mapbz.xc = mapbz.xc - dx0         ;Alignment by comparing the position of the solar limbs
mapbz.yc = mapbz.yc - dy0
mapbx = mapb
mapbx.data = rotate(mapb.data*sin(mapi.data*!dtor)*cos((mapa.data + 270.0)*!dtor),2)  ;The azimuth angle is measured from the CCD+y direction, which is the south, since the solar P angle is ~180 degree
mapbx.roll_angle = 0.0
index = where(mapbx.data le -10000.0)
mapbx.data[index] = 0.0
mapbx.xc = mapbz.xc
mapbx.yc = mapbz.yc
mapby = mapb
mapby.data = rotate(mapb.data*sin(mapi.data*!dtor)*sin((mapa.data + 270.0)*!dtor),2)  ;The azimuth angle is measured from the CCD+y direction, which is the south, since the solar P angle is ~180 degree
mapby.roll_angle = 0.0
mapby.data[index] = 0.0
mapby.xc = mapbz.xc
mapby.yc = mapbz.yc
xc = mapbz.xc
yc = mapbz.yc

intimei = utc2int(mapbz.time)
mjdi  = intimei.mjd + fix((intimei.time + 120000.0)/(24.0*3600.0*1000.0))
timei = (intimei.time + 120000.0) mod (24.0*3600.0*1000.0)
intimei.mjd = mjdi
intimei.time= timei
utctime = anytim2utc(intimei,/vms)
mapbz.time = utctime
mapbx.time = utctime
mapby.time = utctime

loadct,0
wdef,1,800,800
!p.background = 255
plot_map,mapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
write_png,filename+'\'+'bz.png',tvrd()
wdef,2,800,800
plot_map,mapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
plot_vmap,/over,mapbx,mapby,mapbz=mapbz,limit=50,scale=0.1,iskip=30,jskip=30,$
  v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/No_arrow_head  ;,/Noaxis
write_png,filename+'\'+'bxyz.png',tvrd()
print,size(mapbz.data)

sub_map,mapbz,smapbz,xrange=[-75.0,525.0],yrange=[-400.0,-100.0]
print,'Flux balance coefficient:',total(smapbz.data)/total(abs(smapbz.data))
sub_map,mapbx,smapbx,ref_map=smapbz
sub_map,mapby,smapby,ref_map=smapbz
wdef,1,800,800
plot_map,smapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
write_png,filename+'\'+'sbz.png',tvrd()
wdef,2,800,800
plot_map,smapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=180,scale=0.012,iskip=15,jskip=15,$
  v_color=255,axis_color=0,/Nolabels,v_thick=2.0 ;,/No_arrow_head  ;,/Noaxis
write_png,filename+'\'+'sbxyz.png',tvrd()
print,size(smapbz.data)


save,smapbx,smapby,smapbz,filename=filename+'\'+'bxyz_submap.sav'

END