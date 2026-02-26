function Arec = IFWTmatrix_m_VP(Adec,lrow, lcol, theta,fsca,fwav)
%----------------------------------------
% AIM: Reconstruct the approximation at the original resolution level (size) 
% from  the scaling and wavelets coeff. obtained by the decompositions

% INPUT:
% fsca,fwav=  normalization factors for the 1step transforms
% theta = localization parameter for the definition of wavelet 
% lrow= number of decompositions operated on the rows
% lcol= number of decompositions operated on the columns
% Adec = cell array of the transformed matrices
% For k=1:kmax --> Adec{k}= D matrix of all scal. and wave. coeff. 
%                              coming out from the k-th decomposition
% ABOUT THE LENGTH:
% By construction size(D) is of the form [3*n,3*m] 
% scaling coefficients = D(1:n,1:m)
% wavelet coefficients = the remaing part of D

% OUTPUT: 
% Arec= reconstructed matrix 
% It has the same dimension of Adec{1} (multiple of 3)

% NB 
% The size of Adec{k} MUST always be always of the kind [3*n, 3*m].
% At the same step k, if the size of the recontructed scaling matrix A is
% not [n,m] then A is cutted (by 1 or 2 rows or columns) to have that size 
%----------------------------------------
% xx1 =[ 9.988908143703326e-01;-4.162203787241647e-01;5.265210965181560e-02];
%         xx2 =[ 1.001599409402679e+00;-1.439781036554719e-02;-3.843991210102363e-04];
%         fsca = fsca * sqrt(xx2(1) + theta*xx2(2) + theta^2*xx2(3));
%         fwav = fwav * sqrt(xx1(1) + theta*xx1(2) + theta^2*xx1(3));
        
kmax=max(lrow,lcol);
A=Adec{kmax}; % scaling and wavelet coeff. from the last decomposition kmax

%trow=theta; tcol=theta; %parameter theta for row and column decomposition
trow=1/3; tcol=1/3;

 % Reconstruction cicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=kmax:-1:1
%------------------------ 
%INITIAL SETTING
%------------------------ 
B=Adec{k}; [N, M]= size(B); 
%------------------- setting the size of the scaling part
    if k<=lcol
        n=N/3;
    else
        n=N; % we do not decompose the columns
    end
    if k<=lrow
        m=M/3;
    else
        m=M; % we do not decompose the rows
    end
    
%---- scaling coefficients (maybe truncated to be the third part of B)
B(1:n,1:m)=A(1:n,1:m); 
 %------------------------ 
 %  Reconstruction
 %------------------------ 
    if k<=lrow
 %-------- the IFWT (1 step reconstr.) of all rows is taken 
        trow=theta;
        %trow=2*m-1-trow;
        %trow=3*trow;
        for i = 1:N
[At(i,1:M), nminus]= IFWT1step_m_VP(B(i,1:m).',B(i, m+1:M).',trow, fsca,fwav);
        end
    else
        At=B;
    end
    if k<=lcol
  %-------- the IFWT (1 step reconstr.) of all columns is taken
  
        tcol=theta;
        %tcol=2*n-1-tcol;
        %tcol=3*tcol;
        
        for j = 1:M
[At(1:N,j), nminus]= IFWT1step_m_VP(At(1:n,j),At(n+1:N, j),tcol,fsca,fwav);
        end
    end  
 %------------------------ 
 % Final Setting for the next reconstruction
 % %------------------------ 
A = At;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the reconstructed matrix
Arec=A;

end
