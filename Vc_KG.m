function [Vc]=Vc_KG(A,S)
%
% the model is determined by A which is regresson matrix in y=Ax;
% Use spectrum to estimate covariance matrix
% of orthoganal,fitted parameters (b). Further, covariance matrix for 
% orginal paramters (c) is obtained by propagation of error.

% reference 
%
% Doran, Kara J., "Addressing the Problem of Land Motion at Tide Gauges" (2009).
% USF Tampa Graduate Theses and Dissertations.
% https://digitalcommons.usf.edu/etd/1616

N=size(A,1);
n=(0:N-1)';
k=(1:N-1);
K=length(k);
cnk=cos(2*pi*n*k/N); snk=sin(2*pi*n*k/N);

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
%
if mod(N,2)==0
    S1=[S;S(end-1:-1:1)];
else
    S1=[S;S(end:-1:1)];
end
Vb=1.0/N^2*reshape((S1'*reshape(Skij,[K,I*I])),[I,I]);
Vc=diag(pinv(R)*Vb*pinv(R)');

end