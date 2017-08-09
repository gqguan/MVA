function [ T, P ] = MVA_pca( X, opt )
%% process of multivariate analysis
%   principal component decomposition 
%   use the SVD algorithm for opt = 'SVD' to calculate scores and loadings 
%   matrices as large as desired
% refer to 
%   A. Smilde, R. Bro and P. Geladi
%   Multi-way Analysis in Chemistry and Related Fields. 
%   (C) 2004 John Wiley & Sons, Ltd ISBN: 0-471-98691-7
%
%   Chemometrics: Data Analysis for the Laboratory and Chemical Plant.
%   Richard G. Brereton
%   Copyright (C) 2003 John Wiley & Sons, Ltd.
%   ISBNs: 0-471-48977-8 (HB); 0-471-48978-6 (PB)
%
% by Guoqiang GUAN 2017/04/27 lvl-1
%                  2017/05/05 lvl-2 
%
global wrkspace
eps = 1.e-4;
dt = 999;
switch opt
    case 'SVD'
        [U, S, V] = svd(X);
        T = U*S;
        P = V';
    case 'NIPALS'
        % the desired number of PC 
        A = wrkspace.A;
        % initial
        T = zeros(size(X, 1), A);
        P = zeros(A, size(X, 2));
        for a = 1:A
            % the column with greatest sum of squares
            t_guess = MVA_pickcol(X, 'GSS');
            while dt > eps
                p = t_guess'*X/(t_guess'*t_guess);
                p = p/sqrt(p*p');
                t_new = X*p';
                dt = MVA_diff(t_new, t_guess, 'PRESS');
                t_guess = t_new;
            end
            T(:,a) = t_new;
            P(a,:) = p;
            X = X-t_new*p;
            dt = 999;
        end
    otherwise
        warning('Unexpected argument in MVA_pca(..., opt)\n');
end
end