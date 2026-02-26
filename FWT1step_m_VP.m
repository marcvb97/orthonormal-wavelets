function [anew, bnew, nplus] = FWT1step_m_VP(a, theta, fsca,fwav)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In case m=0 we set m=1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
% a= vector to transform of length N1 (must be a column vector)
% theta = parameter in ]0,n[ that define the degree parameter m of our wavelets
%         theta <1  --> m=floor(theta*n)
%         theta >=1 --> m= n-floor(theta);
% fsca,fwav=  normalization factors for the 1step transforms (anew*fsca and
%             bnew*fwav are taken as anew and bnew)
%
% OUTPUT:
% anew= scaling coeff. (column vector)
% bnew =wavelet coeff. (column vector)
% nplus -> indicates how many elements must be added to the input vector a
% nplus=0 --> no element is added
% nplus\ne 0 --> N1+nplus is the adjusted lenght of a 
%          the nplus added elements are fixed equal to the last element of a
%
% The length of at=[anew; bnew] satisfies N=3*n
% The length of anew is n
% The length of bnew is 2*n
%------------------
% Initialization
%------------------
%%%%%% Definition on n and N=3n
    N1 = length(a);
    r1=rem(N1,3);
    if r1==0
        N=N1; n=N/3; 
        nplus=0;
    else  % enlarg  a (of 1 or 2 elem.) by replicating the last entry 
        nplus=3-r1;     % number of elements to be added if r1\ne 0
        for k=1:nplus
            a(N1+k)=a(N1);
        end
        N=N1+nplus; n=N/3; 
    end
%%%%%% Definition of m
        if theta<1
            m=floor(theta*n); 
        else
            m=n-floor(theta);
        end
        %%%%%%%%%%%%%%%%%%%%
        %Check critical case e.g. n=9 and theta=0.1 that gives m=0
        if m==0 %| m==1
          m=1;%  m=2;
        end

%------------------
% Decomposition
%------------------
% STEP 1:
    ap = a(2:3:end); % length(ap)=n
    alphaf = dct(ap * pi/N)* sqrt(n/pi);
    
% STEP 2: 
    app = a;   
    app(2:3:end) = 0; %length(app)=N=3n but n are null
    y2 = dct(app * pi/N)* sqrt(N/pi);
 % computing the beta coeff.
    betaf = y2(1:n-m+1);
    for k = n-m+1:n-1
        betaf(k+1) = (m+n-k)/(2*m)*y2(k+1) - (m-n+k)/(2*m)*y2(2*n-k+1);
    end
% STEP 3:
    mu = ones(n,1);
    for h = n-m+1:n-1
        mu(h+1) = (m^2+(n-h)^2)/(2*m^2);
    end
    ab = (alphaf+betaf)./mu;
  % computing scaling coefficients    
    anew =dct(ab,'Type',3)*sqrt(n/pi);
    anew = anew*fsca;
    %   anew = anew*sqrt(3);
% STEP 4:
    c = zeros(N,1);
    c(1:n-m+1) = ab(1:n-m+1);
    for k = n-m+1:n-1
        c(k+1) = ab(k+1)*(m+n-k)/(2*m);
        c(2*n-k+1) = -ab(k+1)*(m-n+k)/(2*m);
    end
    y5 = dct(c,'Type',3)*sqrt(N/pi);
    y5(2:3:end) = 0;
  % Computing the wavelet coefficients
    bnew = app-y5;
    bnew(2:3:end)=[];
   % bnew=bnew*sqrt(3);
   bnew=bnew*fwav;

end

