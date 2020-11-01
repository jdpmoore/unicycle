function [ plo,pup,plevel,clo,cup,clevel,lambda,mu,disp,strain ] = locate_pnts(gra,gamout,xpnt)
% function LOCATE_PNTS find values for points in GAMRA database gamra & gamout
%
% INPUT:
% gamra  - struct for summary info
% gamout - cell arrays for all patches
% xpnt   - locations of points
%        = [ xx; yy; zz ] [3*pntNum]
%
% OUTPUT:
% plo    - patch lower bounds                  [3*pntNum]
% pup    - patch upper bounds                  [3*pntNum]
% plevel - patch level                         [pntNum*1]
% clo    - cell lower bounds                   [3*pntNum]
% cup    - cell upper bounds                   [3*pntNum]
% clevel - cell level                          [pntNum*1]
% lambda - lambda                              [pntNum*1]
% mu     - mu                                  [pntNum*1]
% disp   - displacements                       [3*pntNum]
% strain - strains                             [9*pntNum]
% note: loop through x first, y next, and z last
%
% first created by Lujia Feng Tue Nov  5 15:04:11 SGT 2013
% last modified by Lujia Feng Tue Nov  5 17:42:17 SGT 2013

% initialize
pntNum =    size(xpnt,2);
plo    =   zeros(3,pntNum);
pup    =   zeros(3,pntNum);
plevel = -1*ones(pntNum,1);
clo    =   zeros(3,pntNum);
cup    =   zeros(3,pntNum);
clevel = -1*ones(pntNum,1);
lambda = -1*ones(pntNum,1);
mu     = -1*ones(pntNum,1);
disp   =   zeros(3,pntNum);
strain =   zeros(9,pntNum);
% data
xlo    = gra.patch_extents.xlo;
xup    = gra.patch_extents.xup;
level  = gra.patch_map.level_number;

% loop through points
for ii=1:pntNum
    % find patch that constains point
    xx      = xpnt(:,ii);
    pInd    = check_inside(xlo,xup,xx);
    inlev   = level(pInd);
    pInd    = pInd(inlev==max(inlev));
    if isempty(pInd)
        error('gamra.locate_pnts ERROR: [%f %f %f] is outside the model domain!',xx);
    end
    if length(pInd)>1
        error('gamra.locate_pnts ERROR: [%f %f %f] belongs to >1 patches at the highest level!',xx);
    end
    plo(:,ii)    = xlo(:,pInd);
    pup(:,ii)    = xup(:,pInd);
    plevel(ii)   = level(pInd);
    % find cell of patch that contains point
    xloCell = gamout{pInd}.xlo;
    xupCell = gamout{pInd}.xup;
    cInd    = check_inside(xloCell,xupCell,xx);
    
    % assign values
    clo(:,ii)    = xloCell(:,cInd);
    cup(:,ii)    = xupCell(:,cInd);
    clevel(ii)   = gamout{pInd}.level(cInd);
    lambda(ii)   = gamout{pInd}.lambda(cInd);
    mu(ii)       = gamout{pInd}.mu(cInd);
    disp(:,ii)   = gamout{pInd}.disp(:,cInd);
    strain(:,ii) = gamout{pInd}.strain(:,cInd);
end


    function [ inInd ] = check_inside(xlo,xup,xx)
        % function CHECK_INSIDE(xlo,xup,xx) check if points reside inside extents
        %
        % INPUT:
        % xlo  - patch lower bounds [3*patchNum]
        % xup  - patch upper bounds [3*patchNum]
        % xx   - locations of points
        %      = [ xx; yy; zz ] [3*1]
        %
        % OUTPUT:
        % inInd - index of patches that contain xx
        %
        % first created by Lujia Feng Tue Nov  5 16:58:50 SGT 2013
        % last modified by Lujia Feng Tue Nov  5 17:19:04 SGT 2013
        
        patchNum = size(xlo,2);
        xx       = xx(:,ones(1,patchNum));
        loInd    = xx>=xlo;
        upInd    = xx<xup;
        inside   = sum(loInd+upInd);
        inInd    = find(inside==6);
        
    end
end