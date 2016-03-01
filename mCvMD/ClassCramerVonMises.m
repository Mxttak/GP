classdef ClassCramerVonMises < handle
% Generic class for computation of the modified Cramer -- von Mises
% distance (mCvMD). For general intro to mCvMD, see Uwe D. Hanebeck and
% Vesa Klumpp, "Localized Cumulative Distributions and a Multivariate
% Generalization of the Cramer-von Mises Distance", MFI 2008
%
% This implementation is porvided as is. For academic purposes only.
%
% AUTHOR: Maxim Dolgov, 01. Mar. 2016
%         maxim.dolgov@kit.edu

properties (SetAccess = protected)
  bmax = -1;
end % function properties

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods
  
  %% ClassCramerVonMises
  function s = ClassCramerVonMises(bmax)
    
    if(~isscalar(bmax) || ~isnumeric(bmax) || bmax<=0)
      error('Input parameter <bmax> must be a positive scalar')
    else
      s.bmax = bmax;
    end    
    
  end % function ClassCramerVonMises
  
end % methods public

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Abstract)
  
  %% ComputeDistance
  ComputeDistance(s,d1,d2)
end % methods abstract

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Access = protected)
  
end % methods protected

end % classdef