function [Adec, lrow, lcol] = FWTmatrix1_m_VP(A, theta, kmax, dmin,fsca,fwav)
%%%%%%%there is the check aboutr the case m=theta*n=0
%----------------------------------------
% AIM: Decompose the matrix A kmax times (or less) while the dimensions 
%(n. of rows and n. of col.) of the scaling part to further decompose 
% is > dmin
%------------------------------------------
% INPUT:
% A = matrix to decompose 
% theta = localization degree of the wavelet we used
% kmax = max number of decompositions we can do on both rows and columns
% dmin= min. dimension of both rows and columns to apply the last decomposition
% (dmin+0,1,2)/3 will be the min. number of scaling elements allowed for each
% row and column
% fsca,fwav=  normalization factors for the 1step transforms 
% OUTPUT:
% lrow= number of decompositions operated on the rows
% nrow= number of decompositions oper ated on the columns
% Adec = cell array --> Adec{k}= D matrix of all scal. and wave. coeff. 
%                              coming out from the k-th decomposition
% ABOUT THE LENGTH:
% By construction size(D) is of the form [3*n,3*m] 
% scaling coefficients = D(1:n,1:m)
% wavelet coefficients = the remaing part of D
%----------------------------------------
% Initialization
%------------------
[N1,M1] = size(A);
dmax=max(N1,M1);

lrow=0; lcol=0; kCONT=0; 

%tcol=theta; trow=theta;
%tcol=3^(floor(log(N1)/log(3)));
%trow=3^(floor(log(M1)/log(3)));

% STARTING of the decomposition cicle
%-------------------------------------

while dmax>dmin && kCONT<kmax
    
 kCONT=kCONT+1;  
%--------- Eventually add rows to get Nj=3*nj
 if N1>dmin
    n1=floor(N1/3); r1=rem(N1,3);
    if r1==0
        nplus=0;
        Nj=N1; nj=n1; 
    else
    % enlarge  A (of 1 or 2 rows) by replicating the last row
        nplus=3-r1;
        %nplusRow(k)=nplus;
        Nj=N1+nplus; nj=Nj/3;
        for k=1:nplus
            A(N1+k,:)=A(N1,:);
        end
    end
 else
     Nj=N1; nj=N1; % in this case we do not decompose the column 
     %               since their length is <= 3
 end
 %--------- Eventually add columns to get Mj=3*mj
 if M1>dmin
    m1=floor(M1/3); r1=rem(M1,3);
    if r1==0
        Mj=M1; mj=m1; 
        nplus=0;
    else
    % enlarg  A (of 1 or 2 columns) by replicating the last column
        nplus=3-r1;
        %nplusCol(k)=nplus;
        Mj=M1+nplus; mj=Mj/3;   
        for k=1:nplus
            A(:,M1+k)=A(:,M1);
        end
    end
 else
     Mj=M1; mj=M1;  % we do not decompose the rows
 end
%------------------
% Now we should have:
% size(A)=[Nj,Mj] with Mj=3*mj and/or Nj=3*nj
%------------------
% Decomposition
%------------------

if Nj>dmin
% the FWT (1step decomp.) of all columns of A is taken
    lcol=lcol+1; 
    tcol=theta;
    %tcol=2*nj-1-tcol;
    %tcol=tcol/3;
    for j = 1:Mj
        [anew, bnew, nplus] = FWT1step_m_VP(A(1:Nj,j), tcol,fsca,fwav);
        at=[anew; bnew]; 
        At1(1:Nj,j) = at;
    end
else
% we do not transform A
    At1=A;        
end
   
if Mj>dmin
% the FWT (1step decomp.) of all rows of At1 is taken
    lrow=lrow+1;
    trow=theta;
    %trow=2*mj-1-trow;
    %trow=trow/3;
    for i = 1:Nj
        [anew, bnew, nplus] = FWT1step_m_VP(At1(i,1:Mj).', trow,fsca,fwav);
        at=[anew; bnew]; 
        At(i,1:Mj) = at.';
    end 
else
% we do not transform At1
    At=At1;
end
%----------------------------------------
% Setting for the nex decomposition step
%----------------------------------------
    A=At(1:nj,1:mj); % matrix of scaling coeff to further decompose
    Adec{kCONT}=At(1:Nj,1:Mj);% matrix of all scaling and wav. coeff. at step k
    N1=nj; M1=mj; dmax=max(N1,M1); % new size of A
%----------------------------------------
end

end