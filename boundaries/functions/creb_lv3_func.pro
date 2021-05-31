pro creb_lv3_func,date1,time1,nx_o,ny_o

time_lst = [date1,'_',time1]
filename1 = 'F:\3D-magnetic-con\3d_myrun\02ambiguity_hmi\05projection\'+strjoin(time_lst)

prep = 1
level = 3

file = find_file(filename1+'\'+'smapbxyz_m.sav')

reducer = 4
nn=n_elements(file)
;nn=1    ;for test
for i=0,nn-1 do begin
  restore,file[i]
  dim=size(smapbz.data,/dim)
  bz=smapbz.data
  print,'Flux Balance Coefficient:',total(bz)/total(abs(bz))
  bx=smapbx.data
  by=smapby.data
  ss=size(bz)
  wdef,1,ss[1],ss[2]
  tvscl,bz
  nx=(size(bz))[1] & ny=(size(bz))[2] & nz=Min([nx,ny]) & nd=Min([nx,ny])/8
  print,'n:',nx,ny,nz,nd
  nx_o = nx/4.0
  ny_o = ny/4.0
  ;  help,bz
;  if nx mod 4 ne 0 then bz = bz[0:-(nx mod 4)-1,0:-1] & bx = bx[0:-(nx mod 4)-1,0:-1] & by = by[0:-(nx mod 4)-1,0:-1]
;  if ny mod 4 ne 0 then bz = bz[0:-1,0:-(ny mod 4)-1] & bx = bx[0:-1,0:-(ny mod 4)-1] & by = by[0:-1,0:-(ny mod 4)-1]
  ;  help,bz
  save,bx,by,bz,filename='F:\3D-magnetic-con\3d_myrun\08others\temp\temp.sav'
  multigrid,file='F:\3D-magnetic-con\3d_myrun\08others\temp\temp.sav',level=level,prep=prep
  filename = 'F:\3D-magnetic-con\3d_myrun\03preprocess\'+strjoin(time_lst)+'_lv3'
;  print,'xc,yc',xc,yc
  spawn,'mkdir '+filename
  spawn,'mv nd.ini '+filename
  spawn,'mv grid.ini '+filename
  spawn,'mv allboundaries.dat '+filename
;  nx=(size(bz))[1] & ny=(size(bz))[2]
;  print,nx/4,ny/4
endfor
spawn,'rm F:\3D-magnetic-con\3d_myrun\08others\temp\temp.sav'

END