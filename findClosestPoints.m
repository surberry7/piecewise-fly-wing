% Surabhi Beriwal
% beriwalsurabhi@gmail.com
% January 2018
    
% input: ppform, points along spline for linear approximation, 
%        number of points selected for each polynomial
%        piece, spline number
% output: parameter values of closest point along curve to
%         linear approximation line segment

function [tValues] = findClosestPoints(ppform, orderedPairs, numberOfPts, spline)

derivativeMatrix = [0 2 0; 0 0 1; 0 0 0];
syms t;
T = [t^2; t; 1];


% tValuesZero = zeros(ppform.pieces * (numberOfPts-1)); % can allocate
                                                        % space for tValues 
                                                        % first
                                                        
% tValues = zeros(ppform.pieces, (numberOfPts-1));                                                  

for piece = 1:1:ppform.pieces % loops through all polynomial pieces of
                              % each spline
    
    % generate explicit formulas for each polynomial piece
    curveCoefs = [ppform.coefs(2*piece-1,:); ppform.coefs(2*piece,:)];
    A = curveCoefs * derivativeMatrix;
    curveSlopes = A * T;
    
    % find closest point on curve to line segment
    % (1) compute slope of each line segment
    % (2) find where slope of line segment and derivative of curve 
    %     are equal
    
    index = numberOfPts * piece + (1 - numberOfPts);
    for j = index:(index+numberOfPts-2)
       dy = orderedPairs(j+1,2)-orderedPairs(j,2);
       dx = orderedPairs(j+1,1)-orderedPairs(j,1);
       lineSlope = dy/dx;
       
       if(j==1)
           tValues = [solve(curveSlopes(2,1)/curveSlopes(1,1)==lineSlope,t)];
       else
           tValues = cat(1,tValues,solve(curveSlopes(2,1)/curveSlopes(1,1)==lineSlope,t));
       end
       
   
    end
    
     
end

