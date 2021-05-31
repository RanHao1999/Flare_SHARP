pro read_boundary_nlfff_func,date1,time1,xc,yc,nx,ny

;========================
;parameters
arcsec2cm = 7.3285e7
xc = xc*arcsec2cm
yc = yc*arcsec2cm
dx = 4.0*0.5d*arcsec2cm
dy=dx
;========================

time_lst = [date1,'_',time1]
filename1 = 'F:\3D-magnetic-con\3d_myrun\03preprocess\'+strjoin(time_lst)
filename2 = 'F:\3D-magnetic-con\3d_myrun\Cartesian\'+strjoin(time_lst)

x = 0.0
bx=fltarr(nx,ny)
by=fltarr(nx,ny)
bz=fltarr(nx,ny)
get_lun,u
openr,u,filename1+'_lv3\'+'allboundaries.dat'
for j=0,ny-1 do begin
  for i=0,nx-1 do begin
    readf,u,x
    bx[i,j]=x
    readf,u,x
    by[i,j]=x
    readf,u,x
    bz[i,j]=x
  endfor
endfor
close,u
free_lun,u
tvscl,bz

spawn,'mkdir '+ filename2 + '\data'
spawn,'mkdir '+ filename2 + '\nlfff_boundary'
spawn,'cp '+ filename1+'_lv3\allboundaries.dat'+' '+filename2+'\nlfff_boundary'
spawn,'cp '+ filename1+'_lv3\grid.ini'+' '+filename2+'\nlfff_boundary'
spawn,'cp '+ filename1+'_lv3\nd.ini'+' '+filename2+'\nlfff_boundary'

sizebz=size(Bz)
nx1=sizebz[1]
nx2=sizebz[2]
filename=filename2+'\nlfff_boundary\nlfff_boundary.dat'
openw,lun,filename,/get_lun
writeu,lun,nx1
writeu,lun,nx2
writeu,lun,double(xc)
writeu,lun,double(yc)
writeu,lun,double(dx)
writeu,lun,double(dy)
writeu,lun,double(Bx)
writeu,lun,double(By)
writeu,lun,double(Bz)
free_lun,lun
print,'Bz range (Gauss):', min(Bz),max(Bz)
;wdef,2,800,800
;tvscl,Bz
print,'Computation domain for nonlinear force-free field:'
nx1=nx1-4
nx2=nx2-4
print,'nx1,nx2',nx1,nx2
print,'xc,yc (cm)',xc,yc
print,'dx,dy (cm)',dx,dy
x1=xc-nx1*dx/2
x2=xc+nx1*dx/2
y1=yc-nx2*dx/2
y2=yc+nx2*dx/2
; output
print,'x,y, and z range (10 Mm):'
print,'        xprobmin1=',strtrim(string(x1*1.e-9),2),'d0'
print,'        xprobmax1=',strtrim(string(x2*1.e-9),2),'d0'
print,'        xprobmin2=',strtrim(string(y1*1.e-9),2),'d0'
print,'        xprobmax2=',strtrim(string(y2*1.e-9),2),'d0'
print,'        xprobmin3=',strtrim(string(0.0*1.e-9+0.1),2),'d0'   ; to lift the domain 1 Mm above 0
print,'        xprobmax3=',strtrim(string((y2-y1)*1.e-9),2),'d0'

openw,lun,filename2+'\parameters.txt',width=256,/get_lun
printf,lun,'Bz range (Gauss):', min(Bz),max(Bz)
;wdef,2,800,800
;tvscl,Bz
printf,lun,'Computation domain for nonlinear force-free field:'
printf,lun,'nx1,nx2',nx1,nx2
printf,lun,'xc,yc (cm)',xc,yc
printf,lun,'dx,dy (cm)',dx,dy
printf,lun,'x,y, and z range (10 Mm):'
printf,lun,'        xprobmin1=',strtrim(string(x1*1.e-9),2),'d0'
printf,lun,'        xprobmax1=',strtrim(string(x2*1.e-9),2),'d0'
printf,lun,'        xprobmin2=',strtrim(string(y1*1.e-9),2),'d0'
printf,lun,'        xprobmax2=',strtrim(string(y2*1.e-9),2),'d0'
printf,lun,'        xprobmin3=',strtrim(string(0.0*1.e-9+0.1),2),'d0'   ; to lift the domain 1 Mm above 0
printf,lun,'        xprobmax3=',strtrim(string((y2-y1)*1.e-9),2),'d0'
free_lun,lun

END