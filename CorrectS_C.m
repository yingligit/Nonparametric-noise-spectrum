function [C,Cr,Ci,Cri]=CorrectS_C(A)
%
% get correction matrix C for residual spectrum from regression model
% the model is determined by A which is regresson matrix in y=Ax;
%
% C:   Correction matrix
% Cr: "real" part correction matrix; 
% Ci: "imaginary" part correction matrix; 
% Cir: correction matrix for expected covariance between the real and
%      imaginary part of fourier transform, because correlation between the
%      two parts could be changed due to the removing of fitted regression
%      models.

% C=Cr+Ci;
% the two part of residual spectrum is different, and may be correlated;
% that is why we derive Cr, Ci, Cir.

% Assume two true spectrum part are equal and is S/2;
%
% Reference: 
% 1) Zhu, Y., Mitchum, G. T., Doran, K. S., Chambers, D. P., & Liang, X. (2021). 
% Distinguishing between regression model fits to global mean sea level reconstructions.
% Journal of Geophysical Research: Oceans, 126, e2021JC017347. https://doi.org/10.1029/2021JC017347
% 2) Doran, Kara J., "Addressing the Problem of Land Motion at Tide Gauges" (2009).
% USF Tampa Graduate Theses and Dissertations.
% https://digitalcommons.usf.edu/etd/1616

N=size(A,1);
n=(0:N-1)';
k=(1:N-1);
K=length(k);
cnk=cos(2*pi*n*k/N); snk=sin(2*pi*n*k/N);

%------ C ---------------
%
[Q,~]=qr(A,0); % orthogonalization
F=Q'*cnk+1i*Q'*snk; F=transpose(F);
% Skij
I=size(Q,2);
Skij=nan(K,I,I);
for ik=1:K
    % Be cautious, transpose of complex matrix is different from that 
    % of real value
    Skij(ik,:,:)=real(transpose(F(ik,:))*conj(F(ik,:)));
end
% Skii
Skii=sum(F.*conj(F),2);
%
C=diag(1.0-2.0/N*Skii)+ ...
   1.0/N^2*(reshape(Skij,[K,I*I])*reshape(Skij,[K,I*I])');
%------ C ---------------

%------ Cr ---------------
Skij_r=nan(K,I,I);
for ik=1:K
    Skij_r(ik,:,:)=real(F(ik,:))'*real(F(ik,:));
end
% Skii
Skii_r=sum(real(F).^2,2);
Cr=diag(0.5-2.0/N*Skii_r)+ ...
   1.0/N^2*(reshape(Skij_r,[K,I*I])*reshape(Skij,[K,I*I])');
%------ Cr ---------------

%------ Ci ---------------
Skij_i=nan(K,I,I);
for ik=1:K
    Skij_i(ik,:,:)=imag(F(ik,:))'*imag(F(ik,:));
end
% Skii
Skii_i=sum(imag(F).^2,2);
Ci=diag(0.5-2.0/N*Skii_i)+ ...
   1.0/N^2*(reshape(Skij_i,[K,I*I])*reshape(Skij,[K,I*I])');
%------ Ci ---------------

%------ Cir ---------------
Skij_ri=nan(K,I,I);
for ik=1:K
    Skij_ri(ik,:,:)=real(F(ik,:))'*imag(F(ik,:));
end
% Skii
Skii_ri=sum(real(F).*imag(F),2);
Cri=diag(-2.0/N*Skii_ri)+ ...
   1.0/N^2*(reshape(Skij_ri,[K,I*I])*reshape(Skij,[K,I*I])');
%------ Cir ---------------

end