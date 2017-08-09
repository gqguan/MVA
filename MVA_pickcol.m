function [ col ] = MVA_pickcol( X, opt )
%% toolkits of multivariate analysis
%   for opt = 'GSS' 
%   pick a column in X with greatest sum of squares 
% refer to 
%   Chemometrics: Data Analysis for the Laboratory and Chemical Plant.
%   Richard G. Brereton
%   Copyright (C) 2003 John Wiley & Sons, Ltd.
%   ISBNs: 0-471-48977-8 (HB); 0-471-48978-6 (PB)
%
% by Guoqiang GUAN 2017/04/27
%
switch opt
    case 'GSS'
        [~, icol] = sort(diag(X'*X), 'descend');
        col = X(:,icol(1));
    otherwise
        warning('Unexpected argument in MVA_pickcol(..., opt)');
end
end

