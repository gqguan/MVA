function [ Z ] = MVA_scaling( X, opt )
%% preprocess of multivariate analysis
%   scaling within the first mode (for opt = 1) 
%   a matrix (X) is scaled such that each 
%   each row is multiplied by a specific scalar
%   or
%   scaling within the second mode (for opt = 2) 
%   each column is multiplied by a certain scalar
%   as in traditional autoscaling
% refer to 
%   A. Smilde, R. Bro and P. Geladi
%   Multi-way Analysis in Chemistry and Related Fields. 
%   (C) 2004 John Wiley & Sons, Ltd ISBN: 0-471-98691-7
%
% by Guoqiang GUAN 2017/04/25
%
switch opt
    case 1
        alt = 2;
    case 2
        alt = 1;
end
W = diag(1./max(X, [], alt));
switch opt
    case 1
        Z = W*X;
    case 2
        Z = X*W;
    otherwise
        warning('Unexpected argument in PreMA_scaling(..., opt)');
end
clear alt W
end