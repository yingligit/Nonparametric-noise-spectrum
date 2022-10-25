function [amin0,Fmin,amin,amax] = ga_min_sort (fcn,amin,amax,nstop,show,Ns,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)

% function [amin,Fmin] = ga_min (fcn,amin,amax,nstop,show,p1,p2,...)
%
% This routine implements a crude version of a genetic algorithm in order
% to minimize a function, F, that is a possibly nonlinear function of a
% parameter vector, a; i.e., F(amin) is minimum.
%
% The "fcn" input is a string giving the name of a m-file that computes
% a vector of F values given a matrix of a vectors. You must also input
% a range that the a vector is within by specifying minimum and maximum
% values in the amin and amax vectors, respectively. The "nstop" input
% determines the stopping criterion, which is that the maximum change
% in any of the a components normalized to the size of amax-amin is
% less than 10^(-nstop); i.e., nstop = 5 means stop when the a vector
% is stable to one part in 10^5. Note that you should scale the problem
% such that the a vector is known to lie in a finite range! The default
% value for nstop is 5, and to use this value input [] for nstop if
% additional parameters are passed to "fcn".
%
% If the "show" input is given, show the results as the iteration proceeds.
% If the "show" is input as [], do not display the intermediate results.
%
% The parameters p1,p2,p3, etc. (up to 10) are passed to "fcn", with the calls
% being made as F=fcn(A,p1,p2,p3,p4, ...) where A is a matrix of a vectors.
% The A matrix is such that each column is an a vector, and the output
% of your function should be a row vector with the same number of
% elements as A has columns; i.e., the number of vectors to check.
%
% Note carefully - if N is the number of elements in a, this routine
% calls F with 10^N guesses at the a vector at a time. In order to
% run efficiently, you must insure that F computes the appropriate
% output within available memory and as efficiently as possible!
%
% The outputs are the best a (amin) and the value of F (Fmin) obtained.

% C is the Correction matrix of one regression model.

amin=amin(:); amax=amax(:); Na=length(amin); %Ns=min(10^Na,1e4);

if isempty(nstop); nstop=4; end
if isempty(show); show=NaN; end

da=amax-amin; %arange=da;
da=max(da,ones(size(da)));    % Allows fixing some parameters.

ck=inf; np=nargin-6;
nn=1;
%for jj=1:1
while ck>1e-3 && nn<30
   % amax and amin should be in desceding order
   for ia=2:Na
       amax(ia)=min(amax(ia-1),amax(ia));
       amin(Na-ia+1)=max(amin(Na-ia+1),amin(Na-ia+2));
   end
   arange=amax-amin;
   %
   a=amin*ones(1,Ns)+rand(Na,Ns).*[arange*ones(1,Ns)];

   %--------------
   % The condition where it requires a be monotonic 
   % if there is no this requirement, comment the line below
   a=sort(a,1,'descend');
   %a1=a(:,1); % first guess of these guesses
   %--------------
   str=['[F]=',fcn,'(a'];
   if np>0
      for n=1:np
          str=[str,',p',int2str(n)];
      end
   end
   str=[str,');'];

   eval(str)
   
   [F,ndex]=sort(F); a=a(:,ndex); Fmin=F(1);
   
   %if issp==1
       % -------------------------
       % Make sum of spectrum equal to sum of the periodogram
       % -------------------------
       %alfa=alfa(ndex);
       %a=a(:,[1:floor(0.05*Ns)])+repmat(log(alfa(1:floor(0.05*Ns))),[Na,1]);
   %else
       % -------------------------
       % No requirement of equal sum of spectrum
       % -------------------------
       a=a(:,[1:floor(0.2*Ns)]);
   %end
   
   amin=min(a')'; amax=max(a')'; arange=amax-amin;
      
   cka=arange./da; ndex=cka<[1/10^nstop];
   amin(ndex)=a(ndex,1); amax(ndex)=a(ndex,1); arange(ndex)=0;

   if ~isnan(show); disp([amin';amax']); end

   ck=max(arange./da);
   nn=nn+1;
end

amin0=a(:,1);
