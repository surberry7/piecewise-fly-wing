% Surabhi Beriwal
% beriwalsurabhi@gmail.com
% October 2017

% This file 1) reads cp2fun to get spline data in B-spline form
%           2) converts b-splines to ppform
%           3) plots ppform
%           4) calls linear approximation function 
%           5) supplies an array of approximated values
%

%timeit(piecewiseFlyWing)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Call cp2fun to extract data                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[F,LM,Lab,R,F0,scale] = cp2fun('/Users/surberry/WingMorphology/WingMorphology/Data/Down1F/Gen03','wing1000.cp');
hold on;

for i = 1:9
    bspline = F{i};   
    ppform = fn2fm(bspline,'pp');
    if(i == 9)
        approximateLinear(ppform);
    else
        approximateLinear(ppform);
    end
end
