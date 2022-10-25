function [F]=MLEpchip_merit_biv_weight(yknot,xknot,logf,sxr,sxi,Cr,Ci,Cri,wk)
%
% maximum likelihood estimation in log space 
% using monotonic cubic spline interpolation
%
% INPUT:
%
% xknot:  log space knot
% yknot:  log space knot values, size(yknot) = [K,Ns], 
%     where Ns is the experiment number
%     where K is knot number
% logf:  all log frequencies, size(logf)=[floor(N/2),1]
% logsp: log spectrum to be fitted, size(logsp)=[N-1,1]
%        low-high-low, as from fft 
% C : correction matrix of a regression model
%
% OUTPUT:
%
% Fmin: merit function of the maximum likelihood
%

Ns=size(yknot,2);
N=size(Cr,2)+1;

yy=pchip(xknot,yknot',logf');% yy size [Ns,length(logf)]

% extend the modeled spectrum to length N-1 (sysmmetric)
if mod(N,2)==0
    yy=[yy,yy(:,end-1:-1:1)];
else
    yy=[yy,yy(:,end:-1:1)];
end

yy=exp(yy');
%yy=(yy');
%
CrS=Cr*yy; CiS=Ci*yy; CriS=Cri*yy;
wk=repmat(wk,[1,Ns]);
F=nansum(wk.* ...
          ((repmat(sxr.^2,[1,Ns]).*(CiS)+...
          repmat(sxi.^2,[1,Ns]).*(CrS)-...
          2*repmat(sxr.*sxi,[1,Ns]).*(CriS))./ ...
          ((CrS).*(CiS)-(CriS).^2)+ ...
          log((CrS).*(CiS)-(CriS).^2)));

end
