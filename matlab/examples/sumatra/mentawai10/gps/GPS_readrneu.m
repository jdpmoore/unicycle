function [ rneu ] = GPS_readrneu(fileName,scaleIn,scaleOut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GPS_readrneu.m	                                %
% read SITE.rneu format GPS time series files                                   %
% SITE.rneu is my reduced format inherited from usgs raw .rneu                  %
% only reads lines started with numbers like 20101021 or 20101021000000         %
% and ignores others                                                            %
%-------------------------------------------------------------------------------%
% The format of *.rneu is                                                       %
% 1        2         3     4    5    6     7     8                              %
% YEARMODY YEAR.DCML NORTH EAST VERT N_err E_err V_err                          %
%-------------------------------------------------------------------------------%
% The format of *.tdpneu is                                                     %
% 1              2         3     4    5    6     7     8                        %
% YEARMODYHHMMSS J2000.SEC NORTH EAST VERT N_err E_err V_err                    %
% Note: time is GPS time                                                        %
%-------------------------------------------------------------------------------%
% The format of *.sneu is                                                       %
% 1              2      3     4    5    6     7     8                           %
% YEARMODYHHMMSS SECOND NORTH EAST VERT N_err E_err V_err                       %
% Note: time is UTC time & earthquake time is 0                                 %
%       s of sneu means second                                                  %
%-------------------------------------------------------------------------------%
%                                                                               %
% INPUT:                                                                        %
% fileName - *.rneu, *.tdpneu, *.sneu files                                     %
% scaleIn  - scale to meter in input fileName                                   %
%    if scaleIn = [], an optimal unit will be determined                        %
% scaleOut - scale to meter in output rneu                                      %
%    if scaleOut = [], units won't be changed                                   %
%                                                                               %
% OUTPUT: an nx8 array (n is number of data points)                             %
%                                                                               %
% first created based on GPS_readrneu.m by Lujia Feng Tue Aug 30 SGT 2011       %
% added zero errors case lfeng Mon Jul  2 17:43:36 SGT 2012                     %
% added scaleIn & scaleOut lfeng Tue Oct 23 17:20:58 SGT 2012                   %
% modified to read tdpneu or sneu files Wed Jul  6 14:47:14 SGT 2016            %
% last modified by Lujia Feng Wed Jul  6 15:15:14 SGT 2016                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% check if this file exists %%%%%%%%%%
if ~exist(fileName,'file'), error('GPS_readrneu ERROR: %s does not exist!',fileName); end

%%%%%%%%%% check extension %%%%%%%%%%
[ ~,~,ext ] = fileparts(fileName);
if strcmpi(ext,'.rneu')
    stdfst = '^[0-9]{8}';  % fst column of data always has format like 20101021
elseif strcmpi(ext,'.tdpneu') || strcmpi(ext,'.sneu')
    stdfst = '^[0-9]{14}'; % fst column of data always has format like 20101021000000
else
    error('GPS_readrneu ERROR: extension %s must be rneu, tdpneu, or sneu!',ext(2:end));
end

%%%%%%%%%% open the file %%%%%%%%%%
fin = fopen(fileName,'r');

rneu = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start reading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)   
    % read in one line
    tline = fgetl(fin);
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % only read lines started with YEARMODY
    isdata = regexp(tline,stdfst,'match');
    if ~isempty(isdata)
       dataCell = regexp(tline,'\s+','split');
       dataStr  = char(dataCell);
       data = str2num(dataStr(1:8,:));
       rneu = [ rneu;data' ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check unit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if scale in the input file is not specified
if ~isempty(scaleOut)
    if isempty(scaleIn)
        errAvg = mean([ rneu(:,6);rneu(:,7);rneu(:,8) ]);
        valMax = max([ rneu(:,3);rneu(:,4);rneu(:,5) ]); 
        if errAvg>1e-5 && errAvg<0.05  % indicate it's m && zero errors
            scaleIn = 1;
        else
            scaleIn = 1e-3;  % indicate it's mm
        end
    end
    scale = scaleIn/scaleOut;
    rneu(:,3:8) = rneu(:,3:8)*scale;
end

fclose(fin);
