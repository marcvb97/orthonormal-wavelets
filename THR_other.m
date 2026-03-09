function [Icomp, zero_el] = THR_other(Idec, Sdec, keep, flag)
% Compression
%-------------
num=numel(Idec);
if flag==0
    stot=sort(abs(Idec),'descend');
    ii=floor(num*keep);  thr=stot(ii);
    zero_el=numel(find(Idec.*(abs(Idec)<thr)));
    Icomp=Idec.*(abs(Idec)>=thr);
else
    Nsca=Sdec(1,1)*Sdec(1,2); % number of scaling coeff.
    Nwav=num-Nsca; % number of details

    X=Idec(Nsca+1:end); %  vector of all details
    [swav,I]=sort(abs(X),'descend');  % vector of all details in decreasing order

    num*keep-Nsca;

    ii=floor(num*keep-Nsca);
    if ii>Nwav || ii<1
        thr=10^3;
        zero_el=0;
        Icomp=Idec;
    else
        X(I(ii+1:end)) = 0.0;
        zero_el = length(I)-ii;
        % thr=swav(ii);
        % I = find(abs(X)<thr);
        % zero_el = length(I);
        % X(I) = 0.0;
        % zero_el=numel(find(X.*(abs(X)<=thr)));
        % Xcomp=X.*(abs(X)>thr);
        Icomp=[Idec(1:Nsca), X];
    end
end