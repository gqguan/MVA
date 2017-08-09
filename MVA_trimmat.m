function [ trimX ] = MVA_trimmat( X, N, opt )
%% process of multivariate analysis
%   trim the input matrix for given row or column (A)
%   for opt = 'col'
%   Eliminate the given number of column
%   for opt = 'row'
%   Eliminate the given number of row
% refer to 
%   Chemometrics: Data Analysis for the Laboratory and Chemical Plant.
%   Richard G. Brereton
%   Copyright (C) 2003 John Wiley & Sons, Ltd.
%   ISBNs: 0-471-48977-8 (HB); 0-471-48978-6 (PB)
%
% by Guoqiang GUAN 2017/05/06 lvl-2
%
[I, J] = size(X);
switch opt
    case 'col'
        trimX = zeros(I, J-1);      
    case 'row'
        trimX = zeros(I-1, J);
        for i = 1:I
            if (i ~= 1) || (i ~= I)
                trimX(1:N-1,:) = X(1:N-1,:);
                trimX(N:I-1,:) = X(N+1:I,:);
            elseif (i == 1)
                trimX = X(2:I,:);
            elseif (i == I)
                trimX = X(1:I-1,:);
            end
        end
    otherwise
        warning('Unexpected argument in MVA_trimmat(..., opt)\n');
end
end

