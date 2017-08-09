function [ A ] = MVA_pcn( X, opt )
%% process of multivariate analysis
%   determine the principal component number (A) for the given X
%   for opt = 'cp-95'
%   A simple rule might be to reject PCs whose cumulative eigenvalues
%   account for less that a certain percentage (e.g. 5 %) of the data
% refer to 
%   Chemometrics: Data Analysis for the Laboratory and Chemical Plant.
%   Richard G. Brereton
%   Copyright (C) 2003 John Wiley & Sons, Ltd.
%   ISBNs: 0-471-48977-8 (HB); 0-471-48978-6 (PB)
%
% by Guoqiang GUAN 2017/05/05 lvl-2
%                  2017/07/17 add case 'cp-99'
%
global wrkspace
% initial PC number
wrkspace = struct('A', rank(X'*X));
[I, J] = size(X);
switch opt
    case 'rankX'
        A = rank(X'*X);
    case 'cp-95'
        [~, V] = MVA_diff(X, X, 'RSS');
        csV = cumsum(V);
        icsV = find(csV > 0.95);
        A = icsV(1);
        fprintf('Using the criteria of cumulative eigenvalues < 5 percent, retained PCs number: %d\n', A);
    case 'cp-99'
        [~, V] = MVA_diff(X, X, 'RSS');
        csV = cumsum(V);
        icsV = find(csV > 0.99);
        A = icsV(1);
        fprintf('Using the criteria of cumulative eigenvalues < 1 percent, retained PCs number: %d\n', A);       
    case 'cv'      
        PRESS = zeros(rank(X'*X), 1);
        RSS = zeros(rank(X'*X), 1);
        for a = 1:rank(X'*X)
            wrkspace.A = a;
            for i = 1:I
                Xest = zeros(size(X));
                Xsub = MVA_trimmat(X, i, 'row');
%                 [~, P] = MVA_pca(Xsub, 'NIPALS');
                P = pca(Xsub, 'Centered', false);
                ti = X(i,:)*P;
                Xest(i,:) = ti*P';
            end
            PRESS(a) = MVA_diff(Xest, X, 'PRESS');
            RSS(a) = MVA_diff(Xest, X, 'RSS');          
            if a > 1
                fprintf('%4d %8.2e %8.2e %8.4f\n', a, PRESS(a), RSS(a), PRESS(a)/RSS(a-1));
                if PRESS(a)/RSS(a-1) > 1
                    break
                end
            else
                fprintf('%4d %8.2e %8.2e\n', a, PRESS(a), RSS(a));
            end
        end
    otherwise
        warning('Unexpected argument in MVA_pcn(..., opt)\n');
end
end

