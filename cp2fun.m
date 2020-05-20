% Surabhi Beriwal
% beriwalsurabhi@gmail.com
% May 2017
% adapted from David Houle's lab

function [F,LM,Lab,R,F0,scale] = cp2fun(pathname,filename)

% close all; %clear all previously generated figures

% cp2fun reads one CP file (format version 1.0 or 2.0) and extracts label
% and control point information, which it uses to put together B-spline
% functions using Matlab Spline Toolbox' SPMAK.
%
% Inputs are path location and file name of the CP file to be read.
%
% Outputs:
    % F: B-spline functions for wing after rotating it to standardized
        % orientation;
    % LM: 12 landmarks, defined by intersection among vein splines;
    % Lab: Alphanumeric labels for wing, read from CP file tag and formatted
        % as a cell array;
    % R: rotation matrix applied to original control points to get the
        % orientation in F;
    % F0: B-spline functions in original orientation;vv
    % scale: scale factor to convert wing coordinates from pixels to mm. Read
        % straight from CP File tag.

% cv = (integer) index for which vein
% nocurves = (integer) total number of veins
% cp = (cell array) control points from .cp file, separated by which vein
%   originally in... pixels? mm? then gets converted to different units and possible rotated?
% scp0 = (array) all control points in one matrix
% LM = array of landmarks / vein intersections NOTE: hardcoded to ignore extra veins
% cpn = array of number of control points for each vein

if filename
    
    if strcmpi(pathname(end),'/')
        
        pathname = pathname(1:end-1);

    end

    fid = fopen([pathname,'/',filename],'r');
    fseek(fid,0,-1);
    cph = fgetl(fid);

    if strcmpi(cph(1:8),'# Image:') %this is a CP ver 1 header

        cphdr = textscan(cph,'# Image: %s %n %n %n %n %s %s %s %s %s %n Resolution: %n %n\n\n');
        fname = lower(cphdr{1}{1});
        cphdr = [{strtok(filename,'.')},cphdr];

        for c = 7:11

            cphdr{c} = cphdr{c}{1};

        end

        cpheader = cphdr([1:10,12,11]);
        nocurves = cell2mat(textscan(fid,'%f\n\n'));

    else %if strcmpi(cph(1:16),'# CPversion: 2.0') %this is a CP ver 2+ header

        cpheader = cell(1,12);
        cpheader{1} = strtok(filename,'.');
        ver2hdrs = {...

            '# Image:',2;
            '# UserCo',3:6;
            '# Imaged',7;
            '# Date: ',8;
            '# Time: ',9;
            '# Strain',10;
            '# Sex: ',12; %notice this is a 7-element field
            '# ScaleF',11;

%             '# Spline',0;
%             '# IsEdit',0;
%             '# IsFini',0;
%             '# IsSpli',0;
%             '# Spline',0;
%             '# OpenRa',0;
%             '# Struct',0;
%             '# Thresh',0;
%             '# FillSi',0;
%             '# ShortC',0;
%             '# LongCl',0;
%             '# Procru',0;

            };

        stop = 0;

        while ~stop % ~ indicates not

            cph = fgetl(fid);
            if strcmpi(cph(1),'#')

                thishdr = cph(1:8);

                switch thishdr(1:7)

                    case '# Sex: '

                        cpheader{ver2hdrs{strcmpi(ver2hdrs,'# Sex: '),2}} = upper(cph(strfind(cph,': ')+2));

                    case '# UserC'

                        cpheader(ver2hdrs{strcmpi(ver2hdrs,'# UserCo'),2}) = textscan(cph(strfind(cph,': ')+2:end),'%n %n %n %n');

                    case '# Scale'

                        cpheader{ver2hdrs{strcmpi(ver2hdrs,'# ScaleF'),2}} = str2double(cph(strfind(cph,': ')+2:end))*2; %<=== NOTE: 2x factor is because even though the scale of the image and anchor points in V2 are scaled twice as large as the "small" version in V1, the landmarks seem to use the same scale as before. This could be a BUG, and need modifying later

                    otherwise

                        if ~isempty(find(strcmpi(cph(1:8),ver2hdrs(:,1))))

                            cpheader{ver2hdrs{strcmpi(cph(1:8),ver2hdrs),2}} = cph(strfind(cph,': ')+2:end);

                        end

                end

            else

                stop = 1;
                nocurves = str2double(cph); %nocurves = number of curves
                
            end

        end

        fname = lower(cpheader{2});
        cpheader{2} = [pathname,'\',fname];

     end


    if and(~isempty(nocurves),~isnan(nocurves))

        for cv = 1:nocurves

            nopoints = cell2mat(textscan(fid,'%f\n\n')); 
            cp{cv,1} = fscanf(fid,'%f %f',[2 nopoints])'; %first instance of cp

        end
       
    end

    fclose(fid);
    
end

if or(isempty(nocurves),isnan(nocurves)) 

    disp(['File ',filename,' contains no data; specimen was likely REJECTED during splining']);
    F = {}; LM={}; X={}; R=[]; F0={}; scale=[];
    
else

    F = cell(nocurves,1); %function structures
    F0 = cell(nocurves,1);
    LM = zeros(12,2);
    Lab = cpheader;

    %scale to mm

    scale = Lab{11};

    %rotating control points to principal axis (improves visualization), original rotation preserved in separate array

    scp0 = cell2mat(cp); %puts all the control points into one matrix - scp0
    cpn = zeros(length(cp),1); %initialized to zero
    cp0 = cp;

    for cv = 1:length(cp) %length(cp) = 9 so cv goes through each of the veins %cv = index for each vein
      
        cpn(cv) = size(cp{cv},1);

    end

    R = pcacov(cov(scp0)); %PCA on a covariance matrix

    if R(1)<0

        R = R*([-1 0;0 1]);

    end

    scp = scp0*R;
    scp = scp - repmat(mean(scp,1),size(scp,1),1); %centering %repmat repeats copies of array

    for cv = 1:length(cp)

        cp{cv} = scp(1:cpn(cv),:)*scale;
        scp(1:cpn(cv),:) = [];

    end

    if cp{1}(1,2) < cp{1}(end,2) %reflect about the y axis

        flipy = 1;

    else

        flipy = -1;

    end

    if cp{1}(1,1) < cp{1}(ceil(length(cp{1})/2),1) %reflect about the X axis

        flipx = -1;

    else

        flipx = 1;

    end

    % forming b-splines

    for cv = 1:length(cp)

        cp{cv}(:,2) = flipy*cp{cv}(:,2);
        cp{cv}(:,1) = flipx*cp{cv}(:,1);

        F{cv,1} = spmak(aptknt(sort(cp{cv}(:,1)),3),cp{cv}');
        F0{cv,1} = spmak(aptknt(sort(cp0{cv}(:,1)),3),cp0{cv}');

    end

    %landmarks (spline intersections):

    LM = [cp{2}(1,:);...
          cp{3}(1,:);...
          cp{4}(1,:);...
          cp{5}(1,:);...
          cp{6}(1,:);...
          cp{6}(end,:);...
          cp{7}(1,:);...
          cp{7}(end,:);...
          cp{8}(1,:);...
          cp{8}(end,:);...
          cp{5}(end,:);...
          cp{9}(1,:)];

      % Graph splines

      figure;

      hold on; axis equal;

      for cv = 1:length(F)
          
          fnplt(F{cv});
          % fnplt(F{cv},'-y'); % plots all splines in yellow, easier to
          %                             visualize linear approximation

          pp = linspace(min(fnbrk(F{cv},'interval')),max(fnbrk(F{cv},'interval')),size(cp{cv},1));
          xx = fnval(F{cv},pp)'; %projects control points on B-splines for plotting
          % xx = cp{cv}; %plots control points with no projections
          % scatter(xx(:,1),xx(:,2),'*');

      end

      scatter(LM(:,1),LM(:,2),'s');

      hold off;

end