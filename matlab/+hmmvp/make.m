function make()
% function MAKE() compiles and links Key Value File Matlab interface
%
%   hmmvp.kvf
%
% usage:
%   
%   hmmvp.make()
%
% note:
%   
%   for openmp, you may have to edit matlabroot/bin/mexopts.sh and change
%     
%     CC='xcrun  -sdk macosx10.9  clang'
%
%   to
% 
%     CC='/sw/bin/gcc-fsf-4.8'
%
%   and make sure 
%
%     MACOSX_DEPLOYMENT_TARGET='10.9'
%     MW_SDK_TEMP="find `xcode-select -print-path` -name MacOSX10.9.sdk"
%
%   are set to the proper values.
%
%  R2014a started using XML files in place of .sh. 
%  Store mexopts.sh in 
%
%     ~/.matlab/R2014b/
%
%  It will then take precedence over the XML compiler files stored in /Applications/MATLAB_R2014a.app/bin/maci64/mexopts/ 
%
%
%  You will need to update /Applications/MATLAB/matlab and add 
%
%  export PATH="$PATH:/Users/sbarbot/Documents/src/unicycle/hmmvp/bin"
%
%

% Configure options
isunix = 1;
use_omp = 1;

% Not sure I have all the flags right. The UNIX part works on my Ubuntu
% 64-bit 4-CPU system with use_omp = 0 or 1. At least one 32-bit Windows user
% has built the code successfully using this script.

p = './+hmmvp';
flags = '-O -I../hmmvp';

if (use_omp) flags = [flags ' -DUTIL_OMP ']; end
if (isunix)
    if (use_omp)
        flags = [flags [' CFLAGS="\$CFLAGS -fopenmp" CXXFLAGS="\$CXXFLAGS -fopenmp" ',...
            'LDFLAGS="\$LDFLAGS -fopenmp -lgomp "']];
    end
    flags = [flags ' CXXLIBS="\$CXXLIBS -Wall -lmwblas -lmwlapack" '];
    mc = sprintf('mex -outdir %s -largeArrayDims %s', p, flags);
else
    % Definitely not sure about the Windows stuff.
    use_omp = 0;
    flags = [flags ' CXXLIBS="$CXXLIBS -lmwblas -lmwlapack "'];
    if (~isempty(dir(sprintf('%s/mexopts.bat', p))))
        mc = sprintf('mex -f %s/mexopts.bat -outdir %s %s', p, p, flags);
    else
        mc = sprintf('mex -outdir %s %s', p, flags);
    end
end

cmd=sprintf('%s ./+hmmvp/hmmvp.cpp ../hmmvp/src/Hmat.cpp ../hmmvp/src/HmatIo.cpp', mc);
fprintf('%s\n',cmd);
eval(cmd);

end
