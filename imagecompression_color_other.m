% AIM: Image compression using custom wavelets.
% The best result is found by varying theta (see vtheta).
% Depending on the value of flag, all decomposition coefficients 
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).
% initial_transformation indicates if TR2D and ITR2D are used.

close all
clear all

%% Settings
flag = 1;   % 0 = compress all coefficients, 1 = compress details only
threshold_type = 1; % 1 = separately, 2 = combined over the channels

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
% vkeep  = [0.5];
vtheta = [0.0];
% vtheta = [0.5];

wname='db2';
wname='bior3.5';

% normalization factors for the 1step transforms
% NB in decomposing we multiply the output by these factors
% in reconstructing we divide the input by these factors
fsca=sqrt(3); fwav=sqrt(3/2);

%% Load and prepare image
testImage = '../PEXELS300/997704.bmp';
testImage = '../13US/fig01.png'
II = imread(testImage);
imshow(II)
R = double(II(:,:,1));
G = double(II(:,:,2));
B = double(II(:,:,3));

[Nor, Mor] = size(R);

%% Decomposition parameters
Kmax = 9;
dmin = 3;

%% Main loop
RESULTS = zeros((Kmax-4+1) * length(vkeep), 6);
cont = 0;

for kstep = 4:Kmax
    for keep = vkeep
        cont = cont + 1;

        %% First theta: initialise error with first vtheta value
        theta = vtheta(1);

        AR = R;
        AG = G;
        AB = B;
        [IdecR,SdecR]=wavedec2(AR,kstep,wname);
        [IdecG,SdecG]=wavedec2(AG,kstep,wname);
        [IdecB,SdecB]=wavedec2(AB,kstep,wname);
        switch threshold_type
            case 1
                [IcompR, zero_elR] = THR_other(IdecR, SdecR, keep);
                [IcompG, zero_elG] = THR_other(IdecG, SdecG, keep);
                [IcompB, zero_elB] = THR_other(IdecB, SdecB, keep);
                zero_el = zero_elR + zero_elG + zero_elB;
            case 2
                [thr, zero_el, IcompR, IcompG, IcompB] = THR_color(IdecR, IdecG, IdecB,...
                    IscalR, IscalG, IscalB, IwavR, IwavG, IwavB, keep, ...
                    lrow, lcol, flag);
        end
        IrecR = waverec2(IcompR, SdecR, wname);
        IrecG = waverec2(IcompG, SdecG, wname);
        IrecB = waverec2(IcompB, SdecB, wname);

        IoutR = IrecR(1:Nor, 1:Mor);
        IoutG = IrecG(1:Nor, 1:Mor);
        IoutB = IrecB(1:Nor, 1:Mor);
        DR    = abs(double(R) - double(IoutR)).^2;
        DG    = abs(double(G) - double(IoutG)).^2;
        DB    = abs(double(B) - double(IoutB)).^2;
        err  = sum(DR(:)+DG(:)+DB(:)) / (numel(R)+numel(G)+numel(B));

        Rfin      = IoutR;
        Gfin      = IoutG;
        Bfin      = IoutB;
        Ifin      = cat(3, Rfin, Gfin, Bfin);
        theta_fin = theta;
        kmax      = kstep;

        %% Quality metrics
        psnrVAL = 10 * log10(255^2 / err);
        SSIMVAL  = (ssim(R, Rfin)+ssim(G, Gfin)+ssim(B, Bfin))/3;

        num = 3*length([IcompR]);

        RESULTS(cont, :) = [keep, psnrVAL, SSIMVAL, (num - zero_el)*100/(3*Nor*Mor), theta_fin, kmax]

        figure(cont)
        imshow(uint8(Ifin))
        axis off equal
    end
end

%% Display results table
RESULTS_best = zeros(length(vkeep),size(RESULTS,2));
for i = 1:length(vkeep)
    I = find(RESULTS(:,1) == vkeep(i));
    [~,J] = max(RESULTS(I,2));
    RESULTS_best(i,:) = RESULTS(I(J(1)),:);
end

varNames = {'CR', 'PSNR', 'SSIM', '\% retained', 'Best $\theta$', 'Decomp'};

disp('***** OUR CASE *****')
T = array2table(RESULTS_best, 'VariableNames', varNames);
disp(T)
