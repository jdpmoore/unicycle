function [ xp,yp,zp,up ] = extents4patch(xlo,xup,level)
% function EXTENTS4PATCH convert patch extents in 3D from GAMRA to
% cocordinates of 6 faces for patch function in MATLAB.
%
% INPUT:
% xlo   - lower bounds
%       = [ xMin; yMin; zMin ] [3*patchNum]
% xup   - upper bounds
%       = [ xMax; yMax; zMax ] [3*patchNum]
% level - on which level patch resides on [patchNum*1]
% dx    - cell size for levels [levelNum*1]
%
% OUTPUT:
% xp, yp, zp, up - input for patch function in MATLAB
% each column is for one patch
%
% first created by Lujia Feng Fri Nov  1 01:00:50 SGT 2013
% last modified by Lujia Feng Mon Nov  4 12:47:41 SGT 2013

xMin = xlo(1,:); % row vector
yMin = xlo(2,:);
zMin = xlo(3,:);
xMax = xup(1,:);
yMax = xup(2,:);
zMax = xup(3,:);

% bottom face
xpBot = [ xMin; xMax; xMax; xMin ]; % [4*patchNum]
ypBot = [ yMin; yMin; yMax; yMax ]; % [4*patchNum]
zpBot = [ zMin; zMin; zMin; zMin ]; % [4*patchNum]

% top face
xpTop = [ xMin; xMax; xMax; xMin ];
ypTop = [ yMin; yMin; yMax; yMax ];
zpTop = [ zMax; zMax; zMax; zMax ];

% left face
xpLef = [ xMin; xMin; xMin; xMin ];
ypLef = [ yMin; yMax; yMax; yMin ];
zpLef = [ zMin; zMin; zMax; zMax ];

% right face
xpRgt = [ xMax; xMax; xMax; xMax ];
ypRgt = [ yMin; yMax; yMax; yMin ];
zpRgt = [ zMin; zMin; zMax; zMax ];

% front face
xpFrt = [ xMin; xMax; xMax; xMin ];
ypFrt = [ yMin; yMin; yMin; yMin ];
zpFrt = [ zMin; zMin; zMax; zMax ];

% back face
xpBak = [ xMin; xMax; xMax; xMin ];
ypBak = [ yMax; yMax; yMax; yMax ];
zpBak = [ zMin; zMin; zMax; zMax ];

% all faces
xp = [ xpBot; xpTop; xpLef; xpRgt; xpFrt; xpBak ];
yp = [ ypBot; ypTop; ypLef; ypRgt; ypFrt; ypBak ];
zp = [ zpBot; zpTop; zpLef; zpRgt; zpFrt; zpBak ];

% put 6 faces of one patch together
xp = reshape(xp,4,[]);
yp = reshape(yp,4,[]);
zp = reshape(zp,4,[]);

if isempty(level)
    up = ones(size(xp));
else
    level = level'; % convert to row vector
    up    = level(ones(4*6,1),:);
    up    = reshape(up,4,[]);
end
end
