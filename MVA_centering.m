function [ Z ] = MVA_centering( X, opt )
%% preprocess of multivariate analysis
%   centering across the first mode (for opt = 1) 
%   the data (X) are centered by subtracting the column-average 
%   from every element in the column
%   or
%   subtracting the row-average from each element in a row is referred to
%   as centering across the second mode (for opt = 2) 
% refer to 
%   A. Smilde, R. Bro and P. Geladi
%   Multi-way Analysis in Chemistry and Related Fields. 
%   (C) 2004 John Wiley & Sons, Ltd ISBN: 0-471-98691-7
%
% by Guoqiang GUAN 2017/04/25
%
switch opt
    case 1
        vecONE = ones(size(X,1),1);
        Z = X-vecONE*mean(X);
    case 2
        rowONE = ones(1,size(X,2));
        Z = X-mean(X,2)*rowONE;
    otherwise
        warning('Unexpected argument in PreMA_centering(..., opt)');
end
clear vecONE rowONE
end