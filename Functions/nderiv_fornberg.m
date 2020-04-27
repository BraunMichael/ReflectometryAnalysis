function [ du ] = nderiv_fornberg( k, xPts, u)
%NDERIV_FORNBERG Uses the Fornberg method for all points to calculate
%   almost-arbitrarily high-order derivatives. It uses a stencil size of 5
%   points in the interior, and 6 for each of the two boundary points. 
%   Here, Nx = length(xPts), and u(ix) = U(xPts(ix)).
% 
% 
% k, the order of the derivative you wish to compute. 
% xPts, the set of points on which u is evaluated. 
% u, the array of values u(x) which are evaluated on xPts (u(i) := u(xPts(i))).
% 
% Requires/Assumes: 
% * xPts should be sorted in ascending order 
% * k < length(xPts)
% 
% RESTRICTIONS: 
%       k must be less than Nx, but greater than 1
%       length(xPts) must be equal to or greater than 6
% CALCULATING DERIV



Nx = length(xPts);
if k >= Nx
   error('*** length(x) must be larger than k')
end
if Nx~=length(u)
   error('*** length(x) must be equal length(u)')
end


du = zeros(size(u));

% Low boundary 2 points
du(1) = dot( fdcoeffF(k,xPts(1),xPts(1:6)), u(1:6) );
du(2) = dot( fdcoeffF(k,xPts(2),xPts(1:6)), u(1:6) );

% interior points
for ix = 3:Nx-2
   du(ix) = dot( fdcoeffF(k,xPts(ix),xPts(ix-2:ix+2)), u(ix-2:ix+2) );
end

% Upper boundary 2 points
du(Nx-1) = dot( fdcoeffF(k,xPts(Nx-1),xPts(Nx-5:Nx)), u(Nx-5:Nx) );
du(Nx  ) = dot( fdcoeffF(k,xPts(Nx  ),xPts(Nx-5:Nx)), u(Nx-5:Nx) );

end

