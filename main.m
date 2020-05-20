% Surabhi Beriwal
% beriwalsurabhi@gmail.com
% January 2018

% This file 1) calls cp2fun.m to get spline data in B-spline form
%           2) converts b-splines to ppform
%           3) calls approximateLinear.m function 
               % input: ppform of each spline
               % output: points along spline for linear approximation, 
               %         number of points selected for each polynomial
               %         piece
%           4) calls findClosestPoints.m function 
               % input: ppform, points along spline for linear approximation, 
               %        number of points selected for each polynomial
               %        piece, spline number
               % output: parameter values of closest point along curve to
               %         linear approximation line segment
%           5) calls evaluateDistance.m function
               % input: ppform, parameter values of closest point along 
               %        curve to linear approximation line segment, points 
               %        along spline for linear approximation, 
               %        number of points selected for each polynomial
               %        piece
               % output: (INCOMPLETE) should output distance between linear
               %         approximation and curve
%           6) displays runtime of all three functions




% call cp2fun.m function
[F,LM,Lab,R,F0,scale] = cp2fun('/Users/surberry/WingMorphology/WingMorphology/Data/Down1F/Gen03','wing1000.cp');
hold on;

for i = 1:9 % loops through each of the nine b-splines
    disp(i); 
    splineNumber = i;
    bspline = F{i};   
    ppform = fn2fm(bspline,'pp'); 
    
    tic
    [linApproxPoints,numberOfPoints] = approximateLinear(ppform);
    disp(strcat('lin approx times(s): ', num2str(toc)));
    tic
    [tValues] = findClosestPoints(ppform,linApproxPoints, numberOfPoints, splineNumber);
    disp(strcat('find pts   times(s): ', num2str(toc)));
    tic
    evaluateDistance(ppform, tValues, linApproxPoints, numberOfPoints);
    disp(strcat('distance   times(s): ', num2str(toc)));

end