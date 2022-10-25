function [C,Ce,yf]=reg_model(A,y)
%
% y is two dimensinonal variable, [time, space point]
%
% C is the fitting model coefficients
% Ce is the standard deviation of C based on formal error
%
[M,N]=size(y);
C=nan(size(A,2),N); % fitting model coefficients matrix
Ce=nan(size(A,2),N); % fitting model coefficients matrix

for k=1:N
    ok=~isnan(y(:,k))&~isnan(sum(A,2));
    if sum(ok)>1
        C(:,k)=A(ok,:)\y(ok,k);        
        % coefficent error
        invATA=inv(A(ok,:)'*A(ok,:));
        Ce(:,k)=nanstd(y(:,k)-A*C(:,k))*sqrt(diag(invATA));
    end
end
% fitted series
yf=A*C;
end