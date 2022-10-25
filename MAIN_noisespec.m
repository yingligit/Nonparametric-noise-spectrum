function [f,sppchip,spres,b,be0,be]=MAIN_noisespec(t,y,A)
%
% 1) Maximum likelihood estimation of noise spectrum from the residual of
% regression model with our newly developed nonparametric model (monotonic
% cubic spline)
% 2) Estimation of regression coefficient standard deviation based on the estimated
% noise spectrum
% 
% INPUT: t: time 
%        y: time series
%        A: desgin matrix or regressor matrix for the regression model 
%
% OUTPUT: f, frequency 
%         sppchip, estimated noise spectrum
%         spres: periodogram of regresssion residual
%         b, regression coefficent
%         be0, standard deviation of b based on formal error (white noise)
%         be, standard deviation of b based on estimated noise spectrum

%
% Version 1: 06/29/2022
%
% Yingli Zhu 
% email: zhuyingliouc@gmail.com 
%
% How to cite
% Zhu, Y., Mitchum, G. T., Doran, K. S., Chambers, D. P., & Liang, X. (2021). 
% Distinguishing between regression model fits to global mean sea level reconstructions.
% Journal of Geophysical Research: Oceans, 126, e2021JC017347. https://doi.org/10.1029/2021JC017347

N=length(t);
dt=t(2)-t(1);
% frequency, we will consider "dt" in freqeuncy at the end of the code
f=(1:floor(N/2))'/N;
logf=log(f);
Nf=length(logf);

%\-------------------------------/
%  1) regression
%\-------------------------------/

% regression, C: regression coefficent, Ce is standard deviation of C based on formal error
% yf: fitted serie, yres: model residual
[b,be0,yf]=reg_model(A,y); % 
yres=y-yf;
ok=~isnan(yres);
yres=interp1(t(ok),yres(ok),t,'linear'); % fill data gaps with linear interpolation
%\-------------------------------/
% 2) get noise spectrum
%\-------------------------------/

Ns=5e3; % Guess number in genetic algorithm 
K=5; % Spline number of monotonic cubic spline
rho1=linspace(0,1,K);
% adjust the spline knot to include more Fourier freqeuncy numbers
xknot=[rho1(1); rho1(2)+0.125; rho1(3)+0.125/2; ...
       rho1(4)+0.125/4; rho1(5)]*log(N/2)-log(N);  
% merit function weights to give more weight to low frequencies
wk2=log((2:Nf+1)./(1:Nf))';
wk2=wk2/sum(wk2);
%
% error spectrum could be underestimated due to regression,  
% the following matrices take into account the overfitting  
[C,Cr,Ci,Cri]=CorrectS_C(A);
C=C(1:Nf,:); Cr=Cr(1:Nf,:); Ci=Ci(1:Nf,:); Cri=Cri(1:Nf,:);
% FFT
sx=fft(yres);
spres=sx(2:floor(N/2)+1).*conj(sx(2:floor(N/2)+1));
if mod(N,2)==0
    sx(N/2+1)=nan; % if N is even number
end
% real and imaginary part of FFT
% NOTE: we add the negative sign in front of imaginary part because FFT in
% our derivation of noise model has a sign opposite to FFT used in MATLAB 
sxr=real(sx(2:floor(N/2)+1)); sxi=-imag(sx(2:floor(N/2)+1));
delta=(C*ones(N-1,1))./ones(Nf,1); % degree of influence of regression
% The noise model will not work very well at frequencies where regression
% remove most of energy, therefore, we ignore these frequencies
sxr(delta<1e-5)=nan; % if the  
sxi(delta<1e-5)=nan;
%
% correct residual periodogram due to overfitting of regression with 
% iteration method developed by Kara Doran and Gary Mitchum 
% spresc1 and spresc1_sm are used to give lower and upper bound of noise spectrum 
[~,spresc1]=KG_noisespec(A,yres);
spresc1(delta<1e-5)=nan;
spresc1_sm=sp_smooth(spresc1,5);
%
% First guesses of lower and upper bound of noise spectrum at the spline knots 
% Note: noise spectrum is first estimated on log frequency space
minsp=min(log(spresc1_sm))*ones(K,1);
maxsp=max(log(spresc1))*ones(K,1);
[yknotg,~,~,~]=ga_min_sort('MLEpchip_merit_biv_weight',...
    minsp, maxsp,5,[],Ns,xknot,logf,sxr,sxi,Cr,Ci,Cri,wk2);
sppchip=exp(pchip(xknot,yknotg,logf));

%\-------------------------------/
% 3) variance of coefficent based on spectrum and model regression
%\-------------------------------/ 
bevar=Vc_KG(A,sppchip);
be=sqrt(bevar);

% consider "dt" in freqeuncy
f=f/dt; 

end
