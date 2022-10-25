%
% Estimate noise spectrum of global mean sea level based on the residual of
% linear trend model 
%
% More cases can be found in
% Zhu, Y., Mitchum, G. T., Doran, K. S., Chambers, D. P., & Liang, X. (2021). 
% Distinguishing between regression model fits to global mean sea level reconstructions.
% Journal of Geophysical Research: Oceans, 126, e2021JC017347. https://doi.org/10.1029/2021JC017347
%
% OUTPUT: f, frequency 
%         sppchip, estimated noise spectrum
%         spres: periodogram of regresssion residual
%         b, regression coefficent
%         be0, standard deviation of C based on formal error (white noise)
%         be, standard deviation of C based on estimated noise spectrum

load church11_mon.mat gmsl t
ok=t>=1900&t<=2010;
t=t(ok); gmsl=gmsl(ok); %
N=length(t);
A=[ones(N,1) (t-mean(t))]; % regressor matrix of linear trend model

[f,sppchip,spres,b,be0,be]=MAIN_noisespec(t,gmsl,A);

figure
loglog(f,spres,'-b',f,sppchip,'-r','linewidth',2)
legend('Residual periodogram','Noise spectrum')
xlabel('Frequency (cpy)','fontsize',20)
ylabel('Spectra (mm^2/cpy)','fontsize',20)
ylim([1e0,1e10])
set(gca,'fontsize',20)
