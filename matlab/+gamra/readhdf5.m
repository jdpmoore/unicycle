function [ gra,gamout ] = readhdf5(fsumName)
% function READHDF5 read GAMRA output including summary and the associated
% data files
%
% INPUT:
% fsumName - GAMRA HDF5 summary files
% fdatName - GAMRA HDF5 data files
%
% OUTPUT:
% gamra  - struct for summary info
% gamout - cell arrays for all patches
%
% EXAMPLE:
%
%   gamra.readhdf5('gamra/homog/kernel-0089-0000.visit')
%
% first created by Lujia Feng Tue Nov  5 11:43:44 SGT 2013
% last modified by Lujia Feng Tue Nov  5 14:17:35 SGT 2013

[wdir,~,~]=fileparts(fsumName);

% read in the summary info
[ gra ] = gamra.readhdf5_summary(fsumName);

% loop through patches to read in data
gamout    = cell(gra.patchNum,1);
for ii=1:gra.patchNum
    % data
    clevel = gra.patch_map.level_number(ii);
    cpatch = gra.patch_map.patch_number(ii);
    processor = gra.patch_map.processor_number(ii);
    [ lambda,mu,disp,strain ] = gamra.readhdf5_data(wdir,processor,clevel,cpatch);
    gamout{ii}.lambda = lambda;
    gamout{ii}.mu     = mu;
    gamout{ii}.disp   = disp;
    gamout{ii}.strain = strain;
    % extents & center
    xlo    = gra.patch_extents.xlo(:,ii);
    xup    = gra.patch_extents.xup(:,ii);
    [ xloCell,xupCell,levelCell ] = gamra.extents2cell(xlo,xup,clevel,gra.dx);
    [ xcrCell ] = gamra.extents2center(xloCell,xupCell);
    gamout{ii}.cellNum= length(levelCell); % number of cells within one patch
    gamout{ii}.xlo    = xloCell; % cell lower bounds
    gamout{ii}.xup    = xupCell; % cell upper bounds
    gamout{ii}.xcr    = xcrCell; % cell centers
    gamout{ii}.level  = levelCell; % cell levels
end

end