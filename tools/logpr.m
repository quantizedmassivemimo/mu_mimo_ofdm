% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : logpr.m
% -------------------------------------------------------------------------
% Description :
%
%   numerical (more) stable way of computing
%   out = log(normcdf(u)-normcdf(l))
%
%   we use the following tricks:
%
%   1) always operate with large negative numbers, if possible
%   2) user erfc instead of normcdf: normcdf(x)=0.5*erfc(-x/sqrt(2))
%   3) use the scaled error function erfcx if argument x is negative: erfcx(x)*exp(-x^2)=erfc(x)
%   4) pull the smaller term out of the logarithm
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================
  
function out = logpr(l,u)

  % -- sanity check

  if any(l>=u)
    error('lower bin position must be smaller than upper bin position!')
  end
  
  % -- do domain transform
  
  ll = l;
  uu = u;
  
  swapidx = abs(u)>abs(l);  
  ll(swapidx) = -u(swapidx);
  uu(swapidx) = -l(swapidx);
  
  % -- pick right approximation depending on uu
  
  setuu = uu<0;
  out(setuu) = -log(2)-uu(setuu).^2/2+log(erfcx(-uu(setuu)/sqrt(2))) + log( 1-exp(-ll(setuu).^2/2+uu(setuu).^2/2).*erfcx(-ll(setuu)/sqrt(2))./erfcx(-uu(setuu)/sqrt(2)) );
  out(not(setuu)) = log(normcdf(uu(not(setuu)))) + log( 1-normcdf(ll(not(setuu)))./normcdf(uu(not(setuu))) );
   
end
