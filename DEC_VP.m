function [Idec, Iscal, Iwav, lrow, lcol] = DEC_VP(I, theta, kstep, dmin,...
    fsca, fwav)
[N0,M0]=size(I); Iin=I; 
% Eventually enlarge the image to get a multiple of 3
%------------------------------------------------------
% adjust the columns
r0=rem(N0,3); 
 if r0~=0
    % enlarge  A (of 1 or 2 rows) by replicating the last row
        Rnplus=3-r0; %number of rows to be added if r0\ne 0
        for k=1:Rnplus
            Iin(N0+k,:)=Iin(N0,:);
        end
 end
% adjust the columns
 r0=rem(M0,3);
 if r0~=0
    % enlarg  A (of 1 or 2 columns) by replicating the last colum
        Cnplus=3-r0; % number of columns to be added if r0\ne 0
        for k=1:Cnplus
            Iin(:,M0+k)=Iin(:,M0);
        end
 end
 % apply decomposition algortithm
 %--------------------------------------------------
    [Idec,lrow,lcol] = FWTmatrix1_m_VP(Iin,theta,kstep, dmin, fsca,fwav);  
%--------------------------------------------------    
% Setting all wavelet coeff. in a unique vector Iwav and 
%         the scaling part in an unique vector Iscal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kmax=max(lrow,lcol); % real number of decomp. we have done
J=Idec{kmax}; [N,M] = size(J); % matrix resulting from the last decomposition
    if kmax<=lcol && kmax<=lrow
        n=N/3; m=M/3;
        I1=J(1:n, m+1:M);
        I2=J(n+1:N,1:m);
        I3=J(n+1:N,m+1:M);
        Iwav=[I1(:);I2(:);I3(:)];
    elseif kmax<=lcol && kmax>lrow
        n=N/3; m=M;
        I1=J(n+1:N,1:M);
        Iwav=I1(:);
    elseif kmax>lcol && kmax<=lrow
        n=N; m=M/3;
        I1=J(1:n, m+1:M);
        Iwav=I1(:);
    end
    Jscal=J(1:n,1:m); 
    Iscal=Jscal(:); % vector of the approx. coeff. in our case
 
 if kmax>1
  for k=kmax-1:-1:1
        J=Idec{k}; [N,M] = size(J); % matrix resulting from the k-th decomposition
    if k<=lcol && k<=lrow
        n=N/3; m=M/3;
        I1=J(1:n, m+1:M);
        I2=J(n+1:N,1:m);
        I3=J(n+1:N,m+1:M);
        Iwav=[Iwav;I1(:);I2(:);I3(:)];
    elseif k<=lcol && k>lrow
        n=N/3; m=M;
        I1=J(n+1:N,1:M);
        Iwav=[Iwav;I1(:)];
    elseif k>lcol && k<=lrow
        n=N; m=M/3;
        I1=J(1:n, m+1:M);
        Iwav=[Iwav;I1(:)];
    end
  end
 end
end
  