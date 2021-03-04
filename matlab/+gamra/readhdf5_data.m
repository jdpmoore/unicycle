function [ lambda,mu,disp,strain ] = readhdf5_data(wdir,processor,level,patch)
% function READHDF5_DATA read lambda,mu,disp,strain out from GAMRA output
% data files
%
% INPUT:
% fileName - GAMRA data output HDF5 files
% level    - patch level starting from 0
% patch    - patch number starting from 0 for each level
% -----------------------------------------------------------------------------
% FILE_CONTENTS {
%  group      /processor.00000
%  group      /processor.00000/level.00000
%  group      /processor.00000/level.00000/patch.00000
%  dataset    /processor.00000/level.00000/patch.00000/Cell lambda
%  dataset    /processor.00000/level.00000/patch.00000/Cell mu
%  dataset    /processor.00000/level.00000/patch.00000/Displacement.00
%  dataset    /processor.00000/level.00000/patch.00000/Displacement.01
%  dataset    /processor.00000/level.00000/patch.00000/Displacement.02
%  dataset    /processor.00000/level.00000/patch.00000/Fault Correction.00
%  dataset    /processor.00000/level.00000/patch.00000/Fault Correction.01
%  dataset    /processor.00000/level.00000/patch.00000/Fault Correction.02
%  dataset    /processor.00000/level.00000/patch.00000/Strain.00
%  dataset    /processor.00000/level.00000/patch.00000/Strain.01
%  dataset    /processor.00000/level.00000/patch.00000/Strain.02
%  dataset    /processor.00000/level.00000/patch.00000/Strain.03
%  dataset    /processor.00000/level.00000/patch.00000/Strain.04
%  dataset    /processor.00000/level.00000/patch.00000/Strain.05
%  dataset    /processor.00000/level.00000/patch.00000/Strain.06
%  dataset    /processor.00000/level.00000/patch.00000/Strain.07
%  dataset    /processor.00000/level.00000/patch.00000/Strain.08
% }
% h5dump -n processor_cluster.00000.samrai
% -----------------------------------------------------------------------------
%
% OUTPUT:
% lambda                  [cellNum*1]
% mu                      [cellNum*1]
% disp    - displacements [3*cellNum]
% strain  - strains       [9*cellNum]
% note: loop through x first, y next, and z last
%
% first created by Lujia Feng Mon Nov  4 15:06:44 SGT 2013
% last modified by Lujia Feng Tue Nov  5 11:40:16 SGT 2013

% display the file structure
%h5disp(fileName);
%info = h5info(fileName);

% use Matlab high-level access functions
% h5read reads dataset
basename = sprintf('/processor.%05d/level.%05d/patch.%05d/',processor,level,patch);
fileName = sprintf('%s/processor_cluster.%05d.samrai',wdir,processor);

% elastic constants
lambda   = h5read(fileName,[basename 'Cell lambda']); % [cellNum*1]
mu       = h5read(fileName,[basename 'Cell mu']); % [cellNum*1]
% displacements
xdisp    = h5read(fileName,[basename 'Displacement.00']); % [cellNum*1]
ydisp    = h5read(fileName,[basename 'Displacement.01']); % [cellNum*1]
zdisp    = h5read(fileName,[basename 'Displacement.02']); % [cellNum*1]
% convert to [3*cellNum] to be consistent with gamra.readhdf5_summary.m
disp     = [ xdisp'; ydisp'; zdisp' ];
% strains
strain0  = h5read(fileName,[basename 'Strain.00']); % [cellNum*1]
strain1  = h5read(fileName,[basename 'Strain.01']); % [cellNum*1]
strain2  = h5read(fileName,[basename 'Strain.02']); % [cellNum*1]
strain3  = h5read(fileName,[basename 'Strain.03']); % [cellNum*1]
strain4  = h5read(fileName,[basename 'Strain.04']); % [cellNum*1]
strain5  = h5read(fileName,[basename 'Strain.05']); % [cellNum*1]
strain6  = h5read(fileName,[basename 'Strain.06']); % [cellNum*1]
strain7  = h5read(fileName,[basename 'Strain.07']); % [cellNum*1]
strain8  = h5read(fileName,[basename 'Strain.08']); % [cellNum*1]
% convert to [9*cellNum] to be consistent with gamra.readhdf5_summary.m
strain   = [ strain0'; strain1'; strain2'; strain3'; strain4'; strain5'; strain6'; strain7'; strain8' ];
end