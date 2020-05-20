% Surabhi Beriwal
% beriwalsurabhi@gmail.com
% September 2017
    
% input: ppform of each spline
% output: points along spline for linear approximation, 
%         number of points selected for each polynomial
%         piece


function [linApproxPoints,numberOfPts] = approximateLinear(ppform)

axis([-1.5 1 -1 1])

% get coefficient matrix for polynomial curves
matrixOfCp = ppform.coefs;

numberOfPts = 4; % number of points per piece 
linApproxPoints = zeros(ppform.pieces * numberOfPts, 2);
% linApproxPoints = zeros(ppform.pieces, numberOfPts, 2);




% get information about breaks for each polynomial piece
    % NOTE: this construction allows for polynomial pieces to be approximated
    %       only within regions of interest
for i = 1:length(ppform.breaks())
    if(i ~= 1 && i ~= length(ppform.breaks()))
        breaks = cat(2,breaks,ppform.breaks(i));
        breaks = cat(2,breaks,ppform.breaks(i));
    end
    
    if(i == 1)
        breaks = ppform.breaks(i);
    end
    
    if(i == length(ppform.breaks()))
        breaks = cat(2,breaks,ppform.breaks(i));

    end
end

% plot pieces of curve
% (1) values of 't' (the parameter) are chosen based on the number of
% points the user specifies the linear approximation should have 
% (2) x and y coordinates are computed using the matrix of coefficients and
% the values of t
% (3) ordered pairs (x(t), y(t)) are plotted

for i = 1:2:length(matrixOfCp(:,1))
  
    t = linspace(0, breaks(i+1)-breaks(i), numberOfPts);
    x = matrixOfCp(i,1) * t.^2 + matrixOfCp(i,2) * t + matrixOfCp(i,3);
    y = matrixOfCp(i+1,1) * t.^2 + matrixOfCp(i+1,2) * t + matrixOfCp(i+1,3);
    X = cat(2,x);
    % X = 
    Y = cat(2,y);
    orderedPairs = [X; Y]; 
    index = (numberOfPts * i + 2 - numberOfPts) / 2;
    linApproxPoints(index:(index+numberOfPts-1),:) = orderedPairs'; % FIX THIS
    plot(orderedPairs(1,:),orderedPairs(2,:), '-k');
    
end 






