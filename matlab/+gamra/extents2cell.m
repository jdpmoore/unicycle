function [ xloCell,xupCell,levelCell ] = extents2cell(xlo,xup,level,dx)
% function EXTENTS2CELL loop thourgh all the patches and divide patch
% extents in 3D into smaller cells using dx
%
% INPUT:
% xlo   - lower bounds
%       = [ xMin; yMin; zMin ] [3*patchNum]
% xup   - upper bounds
%       = [ xMax; yMax; zMax ] [3*patchNum]
% level - on which level patch resides on [patchNum*1]
% dx    - cell size for levels [3*levelNum]
%
% OUTPUT:
% xloCell, xupCell, levelCell for new cells
%
% first created by Lujia Feng Mon Nov  4 14:31:46 SGT 2013
% last modified by Lujia Feng Tue Nov  5 00:54:03 SGT 2013

xloCell   = [];
xupCell   = [];
levelCell = [];
% loop through patches
for ii=1:length(level)
    % patch extents
    xMin = xlo(1,ii); % row vector
    yMin = xlo(2,ii);
    zMin = xlo(3,ii);
    xMax = xup(1,ii);
    yMax = xup(2,ii);
    zMax = xup(3,ii);
    % patch size
    clevel = level(ii)+1; % c style, so need to add 1
    xSize  = dx(1,clevel);
    ySize  = dx(2,clevel);
    zSize  = dx(3,clevel);
    % cell lower extents
    xlin = xMin:xSize:(xMax-xSize);
    ylin = yMin:ySize:(yMax-ySize);
    zlin = zMin:zSize:(zMax-zSize);
    [ ymesh,xmesh,zmesh ] = meshgrid(ylin,xlin,zlin); % loop through x first, then y, and last z
    ll   = [ reshape(xmesh,1,[]); reshape(ymesh,1,[]); reshape(zmesh,1,[]) ];
    xloCell = [ xloCell  ll ];
    % cell upper extents
    xlin = (xMin+xSize):xSize:xMax;
    ylin = (yMin+ySize):ySize:yMax;
    zlin = (zMin+zSize):zSize:zMax;
    [ ymesh,xmesh,zmesh ] = meshgrid(ylin,xlin,zlin); % loop through x first, then y, and last z
    uu   = [ reshape(xmesh,1,[]); reshape(ymesh,1,[]); reshape(zmesh,1,[]) ];
    xupCell   = [ xupCell  uu ];
    levelCell = [ levelCell; clevel*int32(ones(size(ll,2),1)) ];
end
end