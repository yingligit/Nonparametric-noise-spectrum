% smooth spectrum
% Nsm is the width boxcar smoother
function [Ssm]=sp_smooth(S,Nsm)
%
halfwidth=(Nsm-1)/2;
Ssm=nan(size(S));
K=length(S);
%
for ik=1:K
    idx=ik-halfwidth:ik+halfwidth;
    Ssm(ik)=mean(S(idx(idx>=1&idx<=K)));
end
%
end