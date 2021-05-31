pro make_directory,date1,time1

time_lst = [date1,'_',time1]
filename1 = 'F:\3D-magnetic-con\magnetic_modeling_codes-main_1\example\example_Cartesian_uniform_20150827_0524_amrvac20\'
filename2 = 'F:\3D-magnetic-con\3d_myrun\Cartesian\'+strjoin(time_lst)

spawn,'cp '+filename1+'amrvac.par '+filename2
spawn,'cp '+filename1+'mod_usr.t '+filename2
spawn,'cp '+filename1+'potential\amrvac.par '+filename2+'\potential'
spawn,'cp '+filename1+'potential\mod_usr.t '+filename2+'\potential'
;spawn,'cp '+filename1+'amrvac '+filename2
;spawn,'cp '+filename1+'potential\amrvac '+filename2+'\potential'

END