function [ estC ] = MVA_calib( knownC, knownX, givenX, opt )
%% process of multivariate analysis
%   estimate the C with the given X 
%   use the multiple linear regression for opt = 'MLR'
%   use the principal component regression for opt = 'PCR' 
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
%                  2017/07/17 add cases 'PCR/cp-95/NIPALS'
%                                       'PCR/cp-99/NIPALS'
%
global wrkspace
B = 0;
switch opt
    case 'MLR'
        T = knownX;
        P = 1;
    case 'PCR/SVD'
        [T, P] = MVA_pca(knownX, 'SVD');       
    case 'PCR/NIPALS'
        wrkspace.A = MVA_pcn(knownX, 'rankX');
        [T, P] = MVA_pca(knownX, 'NIPALS');
    case 'PCR/cp-95/NIPALS'
        wrkspace.A = MVA_pcn(knownX, 'cp-95');
        [T, P] = MVA_pca(knownX, 'NIPALS');
    case 'PCR/cp-99/NIPALS'
        wrkspace.A = MVA_pcn(knownX, 'cp-99');
        [T, P] = MVA_pca(knownX, 'NIPALS');
    otherwise
        warning('Unexpected argument in MVA_calib(..., opt)\n');
end
R = pinv(T)*knownC;
S = pinv(R)*P;
B = pinv(S);
estC = givenX*B;
end