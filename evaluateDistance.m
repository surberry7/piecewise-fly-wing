% Surabhi Beriwal
% beriwalsurabhi@gmail.com
% March 2018

function evaluateDistance(ppform, tValues, linApproxPoints, numberOfPts)

% input: ppform, parameter values of closest point along 
%        curve to linear approximation line segment, points 
%        along spline for linear approximation, 
%        number of points selected for each polynomial
%        piece
% output: (INCOMPLETE) should output distance between linear
%         approximation and curve

derivativeMatrix = [0 2 0; 0 0 1; 0 0 0];
syms t;
T = [t^2; t; 1];

for piece = 1:1:ppform.pieces
    
    
    curveCoefs = [ppform.coefs(2*piece-1,:); ppform.coefs(2*piece,:)];
    A = curveCoefs * derivativeMatrix;
    slopes = A * T; 
    
    for j = 1:numberOfPts:length(tValues)
        
        for k = j:(j+numberOfPts-1)
            
            % should it be [curveCoefs(1,1) curveCoefs(2,1) 0]
            point = [vpa(subs(slopes(1,1), t, tValues(j))) vpa(subs(slopes(1,1), t, tValues(j))) 0];
            
            % need to remove derivative matrix here and evaluate point and
            % not slope...
        
            lineStart = linApproxPoints(j+piece-1,:);
            lineStart = cat(2,lineStart,0);
            lineEnd = linApproxPoints(j+piece,:);
            lineEnd = cat(2,lineEnd,0);
            
            a = lineStart - lineEnd;
            b = point - lineEnd;
            distance = norm(cross(a,b)) / norm(a); % should output distance
                                 
        
        end
        
    end
       
end
