function [ gra ] = readhdf5_summary(fileName)
% function READHDF5_SUMMARY read GAMRA output summary files
% level 0 always has only one patch that covers the whole domain
% this level-0 patch always has only one child that covers the whole domain
% level 1 always has only one patch that covers the whole domain too
% Note: index uses C convention!!!
%
% INPUT:
% fileName - GAMRA summary HDF5 files
%
% -----------------------------------------------------------------------------
% h5dump -n summary.samrai
% FILE_CONTENTS {
%  group      /BASIC_INFO
%  dataset    /BASIC_INFO/VDR_version_number
%  dataset    /BASIC_INFO/XLO
%  dataset    /BASIC_INFO/child_array
%  - index in the patch list
%  dataset    /BASIC_INFO/child_array_length
%  dataset    /BASIC_INFO/child_pointer_array
%  - offset, number_parents
%  - offset in child_array
%  - offset = -1 means no child patch
%  dataset    /BASIC_INFO/dx
%  dataset    /BASIC_INFO/grid_type
%  dataset    /BASIC_INFO/material_state_variable
%  dataset    /BASIC_INFO/number_dimensions_of_problem
%  dataset    /BASIC_INFO/number_file_clusters
%  dataset    /BASIC_INFO/number_global_patches
%  dataset    /BASIC_INFO/number_levels
%  dataset    /BASIC_INFO/number_patches_at_level
%  dataset    /BASIC_INFO/number_processors
%  dataset    /BASIC_INFO/number_visit_variables
%  dataset    /BASIC_INFO/parent_array
%  - index in the patch list
%  dataset    /BASIC_INFO/parent_array_length
%  dataset    /BASIC_INFO/parent_pointer_array
%  - offset, number_parents
%  - offset in parent_array
%  - offset = -1 means no parent patch
%  dataset    /BASIC_INFO/patch_names_printf
%  dataset    /BASIC_INFO/ratios_to_coarser_levels
%  dataset    /BASIC_INFO/scaling
%  dataset    /BASIC_INFO/time
%  dataset    /BASIC_INFO/time_of_dump
%  dataset    /BASIC_INFO/time_step_number
%  dataset    /BASIC_INFO/var_cell_centered
%  dataset    /BASIC_INFO/var_names
%  - Displacement, Fault Correction, Cell lambda, Cell mu, Strain
%  dataset    /BASIC_INFO/var_number_components
%  dataset    /BASIC_INFO/var_number_ghosts
%  group      /extents
%  dataset    /extents/Cell lambda-Extents
%  dataset    /extents/Cell mu-Extents
%  dataset    /extents/Displacement.00-Extents
%  dataset    /extents/Displacement.01-Extents
%  dataset    /extents/Displacement.02-Extents
%  dataset    /extents/Fault Correction.00-Extents
%  dataset    /extents/Fault Correction.01-Extents
%  dataset    /extents/Fault Correction.02-Extents
%  dataset    /extents/Strain.00-Extents
%  dataset    /extents/Strain.01-Extents
%  dataset    /extents/Strain.02-Extents
%  dataset    /extents/Strain.03-Extents
%  dataset    /extents/Strain.04-Extents
%  dataset    /extents/Strain.05-Extents
%  dataset    /extents/Strain.06-Extents
%  dataset    /extents/Strain.07-Extents
%  dataset    /extents/Strain.08-Extents
%  dataset    /extents/patch_extents
%  - lower, upper, xlo, xup [dim*patchNum]
%  dataset    /extents/patch_map
%    used to connect with data files
%  - processor_number = 0, file_cluster_number = 0  [patchNum*1]
%  - level_number                                   [patchNum*1]
%  - patch_number: patch number in its level        [patchNum*1]
%}
% -----------------------------------------------------------------------------
%
% OUTPUT:
% gamra - struct
% ------------------------------ Basic Info -----------------------------------
% gamra.dim              scalar
% gamra.varNum           scalar
% gamra.varName          [varNum*1] cell
% gamra.compNum          [varNum*1] vector
% gamra.patchNum         scalar
% gamra.levelNum         scalar
% gamra.levelpatchNum    [levelNum*1] vector
% gamra.ratio            [dim*levelNum] matrix
% gamra.scaling          [varNum*1] vector
% ------------------------------ Parent Topo ----------------------------------
% gamra.parentNum        scalar
% gamra.parent           [parentNum*1]
% gamra.parentPnter      struct (offset,number_parents [patchNum*1])
% ------------------------------ Child Topo -----------------------------------
% gamra.childNum         scalar
% gamra.child            [childNum*1]
% gamra.childPnter       struct (offset,number_children [patchNum*1])
% ------------------------------ Patch ----------------------------------------
% gamra.dx               [dim*levelNum]
% gamra.patch_extents    struct (lower,upper,xlo,xup [dim*patchNum])
% gamra.patch_map        struct (processor_number,file_cluster_number,
%                                level_number,patch_number [patchNum*1])
%
% first created by Lujia Feng Thu Oct 31 09:28:01 SGT 2013
% last modified by Lujia Feng Mon Nov  4 21:30:54 SGT 2013

if ~exist(fileName,'file'), error('gamra.readhdf5_summary ERROR: %s does not exist!',fileName); end
fin = fopen(fileName,'r');

% display the file structure
%h5disp(fileName);
%info = h5info(fileName);

% use Matlab high-level access functions
% h5read reads dataset
% ---------------------- Basic Info ----------------------
gra.dim           = h5read(fileName,'/BASIC_INFO/number_dimensions_of_problem'); % scalar =3
gra.varNum        = h5read(fileName,'/BASIC_INFO/number_visit_variables'); % scalar =5
gra.varName       = h5read(fileName,'/BASIC_INFO/var_names'); % [varNum*1] cell
gra.compNum       = h5read(fileName,'/BASIC_INFO/var_number_components'); % [varNum*1] vector
gra.patchNum      = h5read(fileName,'/BASIC_INFO/number_global_patches'); % scalar
gra.levelNum      = h5read(fileName,'/BASIC_INFO/number_levels'); % scalar
gra.levelpatchNum = h5read(fileName,'/BASIC_INFO/number_patches_at_level'); % [levelNum*1] vector, should sum up to patchNum
gra.ratio         = h5read(fileName,'/BASIC_INFO/ratios_to_coarser_levels'); % [dim*levelNum] matrix, all -1 for the 1st level
gra.scaling       = h5read(fileName,'/BASIC_INFO/scaling'); % [varNum*1] vector

% ---------------------- Parent Topo ----------------------
gra.parentNum     = h5read(fileName,'/BASIC_INFO/parent_array_length'); % scalar
if gra.parentNum>1
    gra.parent     = h5read(fileName,'/BASIC_INFO/parent_array'); % [parentNum*1] vector, index in the patch list
    gra.parentPnter= h5read(fileName,'/BASIC_INFO/parent_pointer_array'); % struct (offset,number_parents [patchNum*1])
else
    gra.parent = []; gra.parentPnter = [];
end

% ---------------------- Child Topo ----------------------
gra.childNum      = h5read(fileName,'/BASIC_INFO/child_array_length'); % scalar
if gra.parentNum>1
    gra.child      = h5read(fileName,'/BASIC_INFO/child_array'); % [childNum*1] vector, index in the patch list
    gra.childPnter = h5read(fileName,'/BASIC_INFO/child_pointer_array');; % struct (offset,number_children [patchNum*1])
else
    gra.child = []; gra.childPnter = [];
end

% ---------------------- Patch ----------------------
gra.dx            = h5read(fileName,'/BASIC_INFO/dx'); % [dim*levelNum] matrix
gra.patch_extents = h5read(fileName,'/extents/patch_extents'); % struct (lower,upper,xlo,xup [dim*patchNum])
gra.patch_map     = h5read(fileName,'/extents/patch_map'); % struct (processor_number,file_cluster_number,level_number,patch_number [patchNum*1] int32)

fclose(fin);

end
