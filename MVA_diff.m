function [ E, V ] = MVA_diff( predictX, knownX, opt )
%% toolkits of multivariate analysis
%   for opt = 'PRESS' 
%   calculate the predicted residual error sum of squares 
%   for opt = 'RMSE'
%   calculate the root mean square error
%   for opt = 'RSS'
%   calculate the residual sum of squares
% refer to 
%   Chemometrics: Data Analysis for the Laboratory and Chemical Plant.
%   Richard G. Brereton
%   Copyright (C) 2003 John Wiley & Sons, Ltd.
%   ISBNs: 0-471-48977-8 (HB); 0-471-48978-6 (PB)
%   for opt = 'NRMSE'
%   calculate the normalized RMSE
% refer to https://en.wikipedia.org/wiki/Root-mean-square_deviation
%
% by Guoqiang GUAN 2017/04/26 lvl-1
%                  2017/05/05 lvl-2
%                  2017/08/08 add opt 'NRMSE'
%
global wrkspace
%   get the principal component number from global variables
A = wrkspace.A;
%   Initial values
E = 0;
V = zeros(A, 1);
[I,J] = size(predictX);
switch opt
    case 'PRESS'
        E = norm(predictX-knownX)^2;
    case 'RMSE'
        for i = 1:I
            for j = 1:J
                E = (knownX(i,j)-predictX(i,j))^2+E;
            end
        end
        E = sqrt(E/I/J);
    case 'NRMSE'
        for i = 1:I
            for j = 1:J
                E = (knownX(i,j)-predictX(i,j))^2+E;
            end
        end
        E = sqrt(E/I/J)/(max(max(knownX))-min(min(knownX)));
    case 'RSS'
        X = knownX;
        [T, ~] = MVA_pca(X, 'NIPALS');
        g = ones(A,1);
        for a = 1:A
            g(a) = T(:,a)'*T(:,a);
        end
        for i = 1:I
            for j = 1:J
                E = X(i,j)^2+E;
            end
        end
        for a = 1:A
            V(a) = g(a)/E;
        end
%         csg = cumsum(g);
        % list the eigenvalues
%         fprintf('%4s %10s %10s %12s %7s\n', 'a', 'g(a)', 'V(a)', 'Cumulative %', 'Remarks');
%         fprintf('---- ---------- ---------- ------------ -------\n');
%         for a = 1:A
%             if (csg(a)/E*100 > 95) && (a ~= 1)
%                 fprintf('%4d %10.4e %10.4e %12.4f %7s\n', a, g(a), V(a), csg(a)/E*100, 'Reject');
%             else
%                 fprintf('%4d %10.4e %10.4e %12.4f %7s\n', a, g(a), V(a), csg(a)/E*100, 'Retain');
%             end
%         end
        E = E-sum(g);
    otherwise
        warning('Unexpected argument in MVA_diff(..., opt)\n');
end
end