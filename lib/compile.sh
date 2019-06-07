# ! /bin/bash


# new addition: 
#rm *.o

rm pygfunc.c

gfortran -fPIC -c pygfunc.f90  
gfortran -fPIC -c write_data_phantom.f90   #phantread_module.f90 

python setup.py build_ext --inplace #>> error.out

#python ../work/run.py ../work/mainsequence.cfg 


#cp star_00000.tmp ../../run_phantom/M2H_test/
#cp  ../../run_phantom/M2H_test/star_00000.tmp ../../run_phantom/M2H_test/star_00000.back