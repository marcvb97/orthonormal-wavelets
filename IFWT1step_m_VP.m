function [aa, nminus] = IFWT1step_m_VP(a,b,theta, fsca,fwav)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In case m=0 we set m=1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT:
% a= scaling coeff. (must be a column vector)
% b =wavelet coeff. of length 2*n (must be a column vector)
% theta = parameter in ]0,n[ that define the degree parameter m of our wavelets
%         theta <1  --> m=floor(theta*n)
%         theta >=1 --> m= n-floor(theta);
% fsca,fwav=  normalization factors for the 1step transforms 
%               (a/fsca and b/fwav are taken as a and b)
% divided by 

% OTPUT: 
% aa= column vector transformed (of length 3*n <= length(a)+length(b))

% NB 
% The length of b MUST always be even (2*n). If the length(a) is not n 
% then a is cutted (by 1 or 2) to have n entries. In this case:
% length(aa)<length(a)+length(b) 
%------------------
% Initialization
%------------------
%%%%%% Definition on n and N=3n
n = length(b)/2; % n=floor(length(b)/2);
nminus=length(a)-n; % number of elements we cut to have
%                       length(a)=2*length(b) 
N=3*n; 
%%%%%% Definition of m
if theta<1
    m=floor(theta*n); 
else
    m=n-floor(theta);
end
%%%%%%%%%%%%%%%%%%%%
        %Check critical case e.g. n=9 and theta=0.1 that gives m=0
        if m==0 %| m==1
           m=1; % m=2;
        end
%%%%%%%%%%%%%%%%%%%%%
anew = a(1:n); % scaling coefficients (maybe truncated)
bnew=b(1:2*n); % wavelets coefficients
% NEW normalization 
anew=anew/fsca;
bnew=bnew/fwav;  
%----------------
% Reconstruction
%----------------

% STEP 1: (computing alpha coeff.)
    alphaf2 = dct(anew * pi/n)* sqrt(n/pi);

% STEP 2:
    b2 = zeros(N,1);
    for k = 1:n
        ib2=1+3*(k-1); ib=1+2*(k-1);
        b2(ib2) = bnew(ib);
        b2(ib2+2) = bnew(ib+1);
    end
    y2=dct(b2* pi/N)* sqrt(N/pi);
  % computing beta coefficients
    betaf2 = y2(1:n-m+1);
    for k = n-m+1:n-1
        betaf2(k+1) = (m+n-k)/(2*m)*y2(k+1) - (m-n+k)/(2*m)*y2(2*n-k+1);
    end
    
% STEP 3: (compute ap2= a-prime coeff.)
    ap2 = dct(3*betaf2,'Type',3)* sqrt(n/pi);
    ap2 = anew - ap2;

% STEP 4: (compute app2= a-douple-prime coeff)
    c = zeros(N,1);
    c(1:n-m+1) = alphaf2(1:n-m+1);
    for k = n-m+1:n-1
        c(k+1) = alphaf2(k+1)*(m+n-k)/(2*m);
        c(2*n-k+1) = -alphaf2(k+1)*(m-n+k)/(2*m);
    end
    y5 = dct(c,'Type',3)* sqrt(N/pi);
    y5(2:3:end) = 0;
    app2 = b2 + y5;

% COMPUTE THE RECONSTRUCTED OUTPUT 
    aa = app2;
    aa(2:3:end) = ap2;
    %aa=aa/sqrt(3);

end
