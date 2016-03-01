classdef ClassDiracDiracDistance < ClassCramerVonMises
% Implementation of the modified Cramer -- von Mises distance (mCvMD)
% between two Dirac densities as proposed in Uwe D. Hanebeck, "Optimal
% Reduction of Multivariate Dirac Mixture Densities",
% at-Automatisierungstechnik, 2015
%
% This implementation is porvided as is. For academic purposes only.
%
% AUTHOR: Maxim Dolgov, 01. Mar. 2016
%         maxim.dolgov@kit.edu
%
% Uses implementation of the Euler gamma function by Paul Godfrey,
% pgodfrey@conexant.com available on Matlab Central.
%
% TODO: - implement input checks
%       - implement method AutoProduct
  
methods
  
  %% ClassDiracDiracDistance
  function s = ClassDiracDiracDistance(bmax)
    s@ClassCramerVonMises(bmax);
  end % function ClassDiracDiracDistance
  
  %% ComputeDistance
  function [dist,means] = ComputeDistance(s,d1,d2)
    [m1,n1] = size(d1);
    [m2,n2] = size(d2);
    
    if(m1 ~= m2)
      error('Input parameters <d1> and <d2> must have compatible dimensions')
    elseif(~isnumeric(d1) || ~isnumeric(d2))
      error('Input parameters <d1> and <d2> must be numeric')
    end
  
    means(:,1) = mean(d1,2);
    means(:,2) = mean(d2,2);
    
    dist = pi^(m1/2)*(s.CrossProduct(d1,d1,n1,n1) - 2*s.CrossProduct(d1,d2,n1,n2) + s.CrossProduct(d2,d2,n2,n2) + 2*(log(4*s.bmax^2)-s.EulerGamma)*(means(:,1)-means(:,2))'*(means(:,1)-means(:,2)) )/8;
  end % function ComputeDistance
  
end % methods public

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Access = private)

%% AutoProduct
  function out = AutoProduct(s,d,n)
  end % function AutoProduct
  
  %% CrossProduct
  function out = CrossProduct(s,d1,d2,n1,n2)
    ind = allcomb(1:n1,1:n2);
    tmp = (d1(:,ind(:,1))-d2(:,ind(:,2)));
    out = sum(s.Xlog(sum(tmp.*tmp,1)))/(n1*n2);
  end % function CrossProduct

end % methods private

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
methods (Static)
    
  %% Eulergamma
  function [eg] = EulerGamma()
  %Euler-Mascheroni constant = -Psi(1) = 0.5772156649015328606...
  %
  %see also: psi, psiz
  %

  %Paul Godfrey
  %pgodfrey@conexant.com
  %8-7-2000

    eg=0.5772156649015328606;

  end % function Eulergamma
  
  %% Xlog
  function y = Xlog(x)

    y = real(x.*log(x));
    y(isnan(y)) = 0;
  end % function Xlog
  
end % methods static 

end % classdef