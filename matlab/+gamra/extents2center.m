function [ xcr ] = extents2center(xlo,xup)
% function EXTENTS2CENTER find the centers given the extents of patches
%
% INPUT:
% xlo   - patch lower bounds
%       = [ xMin; yMin; zMin ] [3*patchNum]
% xup   - patch upper bounds
%       = [ xMax; yMax; zMax ] [3*patchNum]
%
% OUTPUT:
% xcr   - patch center
%       = [ xMid; yMid; zMid ] [3*patchNum]
%
% first created by Lujia Feng Tue Nov  5 10:32:29 SGT 2013
% last modified by Lujia Feng Tue Nov  5 10:37:55 SGT 2013

xcr = xlo + 0.5*(xup - xlo);

end