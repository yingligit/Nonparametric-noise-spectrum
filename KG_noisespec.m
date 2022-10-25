function [Sres,S]=KG_noisespec(A,y)
% get corrected spectrum using the iteration method proposed by Kara Doran and Gary Mitchum  
% with parameter errors
%
% Input:
% A:       The basis function of regression model
% y:       The residual from the regression model
%
% Output:
% Sres:    Residual spectrum 
% S:       Corrected residual Spectrum using the iteration method
%-----
% NOTE:: The frequency for Sres S doesn't include zero frequency 
%-----
% reference 
%
% Doran, Kara J., "Addressing the Problem of Land Motion at Tide Gauges" (2009).
% USF Tampa Graduate Theses and Dissertations.
% https://digitalcommons.usf.edu/etd/1616

N=length(y);
% Fourier transform of othogonalized basis function
n=(0:N-1)';
k=(1:N-1);
K=length(k);
cnk=cos(2*pi*n*k/N); snk=sin(2*pi*n*k/N);

% get power spectrum
Xk=y'*cnk+1i*y'*snk;
Sres=Xk.*conj(Xk);
Sres=Sres(:);
%
[Q,R]=qr(A,0); % orthogonalization
F=Q'*cnk+1i*Q'*snk; F=transpose(F);
% Skij
I=size(Q,2);
Skij=nan(K,I,I);
for ik=1:K
    % Be cautious, transpose of complex matrix is different from that 
    % of real value
    Skij(ik,:,:)=real(transpose(F(ik,:))*conj(F(ik,:)));
end
Skii=sum(F.*conj(F),2);

S=N/2*ones(K,1);
% iterate to get true spectrum
for p=1:2
    gama1=1.0 ...
        +1.0/N^2*(reshape(Skij,[K,I*I])*reshape(Skij,[K,I*I])'*S)./S ...
        -2.0/N.*Skii;
    S=1./gama1.*Sres;
end
    
% output the first half frequency band spectrum with floor(N/2)
% frequencies.
Sres=Sres(1:floor(N/2));
S=S(1:floor(N/2));

end
