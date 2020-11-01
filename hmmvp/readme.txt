# compiling hmmvp and installing the Matlab interface.

1) open matlabroot/bin/mexopts.sh and change

       CC='xcrun -sdk macosx10.X clang'

    to 

       CC='/sw/bin/gcc-fsf-X'

    where 'X' means depend on what version you have, check by ls /sw/bin/gcc-fsf*

2) make sure 

       MACOSX_DEPLOYMENT_TARGET='10.X'
       MW_SDK_TEMP="find `xcode-select -print-path` -name MacOSX10.X.sdk"

    are set to the proper values. The 'X' above means version of your OS, go to Apple-> 'about this Mac' to check the version.


3) open Makefile and make sure 

       CC=/sw/bin/g++  

    is consistent with g++ you have, for me it is g++-5 in '/sw/bin/', so I set the above CC to 
       CC=/sw/bin/g++-5


4) go to /Application/MATLAB/bin/ and open matlab then add
    
       export PATH="$PATH:/Users/rino/Documents/src/unicylce/hmmvp/bin"


5) build the mex files inside Unicycle

      cd PATH_TO_UNICYCLE/matlab
      hmmvp.make()

   update +hmmvp/make.m to fit your environment. 
