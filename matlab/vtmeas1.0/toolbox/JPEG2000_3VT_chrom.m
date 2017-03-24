% 3AFC test codes for chrominance components
% Author: Feng Liu, the University of Arizona, 3/2/2017

fprintf('\nVERY IMPORTANT! PLEASE NOTE BEFORE PROCEEDING:\n\nMeasurements of VTs at lower JND levels has to be totaly finished before VTs at higher JND levels are measured.\n\nVTs for HH1-5, HL1-5, LH1-5 and LL5 are needed at each JND level for the encoder. The VTs\nof HL and LH subbands at a JND level are very close to each other. So the final results for\nthese two subbands are pooled by taking the minimum.\n\nThe VT do not vary a lot as the subband coefficient variance changes. So a single typical\nvariance from [1] is used for each subband.\n\n32 3AFC trials with the question "Which stimulus is different from the other two?" are taken\nfor each subband at a given subband coefficient variance and JND level. 10-second observation\ntime limit is enabled for each trial. The subject has to answer the question by pressing\n"1", "2" or "3" at the end of each trial.\n\nFor more details of the measurement experiment set-up and procedure, as well as how to use\nthe VTs in the encoder, please refer to [2]\n\n[1] H. Oh, A. Bilgin and M. W. Marcellin, "Visually Lossless Encoding for JPEG2000," in\nIEEE Transactions on Image Processing, vol. 22, no. 1, pp. 189-201, January, 2013.\n\n[2] Feng Liu, Yuzhang Lin, Eze L. Ahanonu, Michael W. Marcellin, Amit Ashok, Elizabeth A. Krupinski\nand Ali Bilgin, "Visibility Thresholds for Visually Lossy JPEG2000," Proc. SPIE 9971,\nApplications of Digital Image Processing XXXIX, 99711P, September, 2016.\n\n');
%% test parameters

B=8; % only 8-bit is supported for each chrominance.
VTidx=input('JND Level (Positive Integer): ');
cont=1;
while(cont==1)
    if(isempty(VTidx))
        VTidx=input('JND Level (Positive Integer): ');
    elseif((VTidx<1)||(floor(VTidx)-VTidx~=0))
        VTidx=input('JND Level (Positive Integer): ');
    else
        cont=0;
    end
end
level=input('Wavelet Decomposition Level (1~5): ');
cont=1;
while(cont==1)
    if(isempty(level))
        level=input('Wavelet Decomposition Level (1~5): ');
    elseif((level<1)||(level>5)||(floor(level)-level~=0))
        level=input('Wavelet Decomposition Level (1~5): ');
    else
        cont=0;
    end
end
band=input('Subband Direction (0-LL, 1-HL, 2-LH, 3-HH): ');
cont=1;
while(cont==1)
    if(isempty(band))
        band=input('Subband Direction (0-LL, 1-HL, 2-LH, 3-HH): ');
    elseif((band<0)||(band>3)||(floor(band)-band~=0))
        band=input('Subband Direction (0-LL, 1-HL, 2-LH, 3-HH): ');
    else
        cont=0;
    end
end
loadhead={'LL','HL','LH','HH'};
blk_size_m = 64;   %default codeblock size (m=height, n=width), actual codeblock sizes are internally adjusted according to the size of subbands
blk_size_n = 64;
sig2 = [0.18 0.17;1.33 1.43;0.72 0.74;3.06 3.52;1.16 1.20;4.26 5.14;1.43 1.37;5.34 7.85;1.52 1.45;5.34 8.93;150.08 109.51]; % typical variance values
sig2=sig2((level-1)*2+(3-ceil(band/2)),comp);

%% test setup
if(comp==1)
    if(exist([pwd '/chromthresh/' loadhead{band+1} sprintf('%dCbvt%d_3afc.mat',level,VTidx)])~=0)
        cd ..;
        error('Enough measurements has been made for this subband, no need to measure any more.');
    end
    if(VTidx~=1)
        if(band~=0)
            if(exist([pwd '/data/' sprintf('CbVT%d.mat',VTidx-1)])~=0)
                load([pwd '/data/' sprintf('CbVT%d.mat',VTidx-1)]);
            else
                findCbVT(VTidx-1);
            end
            delta_base=T(3-ceil(band/2),level);
        else
            if(exist([pwd '/data/' sprintf('LL5Cbvt%d_3afc.mat',VTidx-1)])~=0)
                load([pwd '/data/' sprintf('LL5Cbvt%d_3afc.mat',VTidx-1)]);
            else
                findCbVT(VTidx-1);
            end
            delta_base=t;
        end
    else
        delta_base=0;
    end
else
    if(exist([pwd '/chromthresh/' loadhead{band+1} sprintf('%dCrvt%d_3afc.mat',level,VTidx)])~=0)
        cd ..;
        error('Enough measurements has been made for this subband, no need to measure any more.');
    end
    if(VTidx~=1)
        if(band~=0)
            if(exist([pwd '/data/' sprintf('CrVT%d.mat',VTidx-1)])~=0)
                load([pwd '/data/' sprintf('CrVT%d.mat',VTidx-1)]);
            else
                findCrVT(VTidx-1);
            end
            delta_base=T(3-ceil(band/2),level);
        else
            if(exist([pwd '/data/' sprintf('LL5Crvt%d_3afc.mat',VTidx-1)])~=0)
                load([pwd '/data/' sprintf('LL5Crvt%d_3afc.mat',VTidx-1)]);
            else
                findCrVT(VTidx-1);
            end
            delta_base=t;
        end
    else
        delta_base=0;
    end
end
trialsDesired=32;   %QUEST trials
pThreshold=0.82; %QUEST correct-point
%--------------------------------------------------------------

err = struct('xabs_max',0,'yabs_max',0,'xabs_mean',0,'yabs_mean',0,'yabs_big_mean',0,'delta',0);
% Provide the prior knowledge to QuestCreate, and receive the data struct "q".
tmax=128;
tGuess=input(sprintf('Estimate threshold (positive value ranging in 0-%d, e.g., 1.5): ',tmax));
cont=1;
while(cont==1)
    if(isempty(tGuess))
        tGuess=input(sprintf('Estimate threshold (positive value ranging in 0-%d, e.g., 1.5): ',tmax));
    elseif((tGuess<=0)||(tGuess>tmax))
        tGuess=input(sprintf('Estimate threshold (positive value ranging in 0-%d, e.g., 1.5): ',tmax));
    else
        cont=0;
    end
end
tGuessSd=input(sprintf('Estimate the standard deviation of your guess, above (positive value ranging in 0-%d, e.g., 3.0): ',tmax));
cont=1;
while(cont==1)
    if(isempty(tGuessSd))
        tGuessSd=input(sprintf('Estimate the standard deviation of your guess, above (positive value ranging in 0-%d, e.g., 3.0): ',tmax));
    elseif((tGuessSd<=0)||(tGuessSd>tmax))
        tGuessSd=input(sprintf('Estimate the standard deviation of your guess, above (positive value ranging in 0-%d, e.g., 3.0): ',tmax));
    else
        cont=0;
    end
end

beta=3.2;delta=0.01;gamma=0.333;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,0.01,2*tGuessSd);
q.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

% On each trial we ask Quest to recommend an intensity and we call QuestUpdate to save the result in q.

wrongRight={'wrong','right'};

%% VT measurement starts

for k=1:trialsDesired
	% Get recommended level.  Choose your favorite algorithm.
	tTest=QuestQuantile(q);	% Recommended by Pelli (1987), and still our favorite.
% 	tTest=QuestMean(q);		% Recommended by King-Smith et al. (1994)
% 	tTest=QuestMode(q);		% Recommended by Watson & Pelli (1983)
	
    if tTest < delta_base
        tTest = delta_base;
    end
    
    if(VTidx==1)
        [response, err] = JND_trial_3VT(B, comp, level, band, blk_size_m, blk_size_n, 0, tTest, sig2);
    else
        [response, err] = JND_trial_3VT(B, comp, level, band, blk_size_m, blk_size_n, delta_base, tTest, sig2);
    end
    
    tTestlist(k) = tTest;
    responselist(k) = response;
    errlist(k) = err;
        
    
	fprintf('Trial %3d at %4.5f is    %s : xmax-%4.5f, ymax-%4.5f, xavg-%4.5f, yavg-%4.5f, ybavg-%4.5f\n',k,tTestlist(k),char(wrongRight(responselist(k)+1)),err.xabs_max,err.yabs_max,err.xabs_mean, err.yabs_mean, err.yabs_big_mean);
    
	% Update the pdf
	q=QuestUpdate(q,tTest,response); % Add the new datum (actual test intensity and observer response) to the database.
end
    
% Ask Quest for the final estimate of threshold.
t=QuestMean(q);		% Recommended by Pelli (1989) and King-Smith et al. (1994). Still our favorite.
sd=QuestSd(q);
fprintf('Final threshold estimate (mean) is %.4f +- %.4f\n',t,sd);
% t=QuestMode(q);	% Similar and preferable to the maximum likelihood recommended by Watson & Pelli (1983). 
% fprintf('Mode threshold estimate is %4.2f\n',t);
%fprintf('\nYou set the true threshold to %.2f.\n',tActual);
fprintf('Quest knew only your guess: %.2f +- %.2f.\n',tGuess,tGuessSd);

%% data saving and synthesis

if(comp==1)
    save([pwd '/chromthresh/' loadhead{band+1} sprintf('%dCbvt%d_3afc.mat',level,VTidx)],'t','sd'); % Cb Component
    bandsign=zeros(3,5);
    for kz=1:5
        for kj=1:3
            if(exist([pwd '/chromthresh/' loadhead{kj+1} sprintf('%dCbvt%d_3afc.mat',kz,VTidx)])~=0)
                bandsign(kj,kz)=1;
            else
                fprintf(['VT measurement for Cb ' loadhead{kj+1} sprintf('%d subband at JND=%d is still needed.\n',kz,VTidx)]);
            end
        end
    end
    if(exist([pwd '/chromthresh/' sprintf('LL5Cbvt%d_3afc.mat',VTidx)])==0)
        fprintf('VT measurement for Cb LL5 at JND=%d is still needed.\n',VTidx);
    elseif(sum(bandsign(:))==15)
        findCbVT(VTidx);
        extrafit_Cb(VTidx);
    end
else
    save([pwd '/chromthresh/' loadhead{band+1} sprintf('%dCrvt%d_3afc.mat',level,VTidx)],'t','sd'); % Cr Component
    bandsign=zeros(3,5);
    for kz=1:5
        for kj=1:3
            if(exist([pwd '/chromthresh/' loadhead{kj+1} sprintf('%dCrvt%d_3afc.mat',kz,VTidx)])~=0)
                bandsign(kj,kz)=1;
            else
                fprintf(['VT measurement for Cr ' loadhead{kj+1} sprintf('%d subband at JND=%d is still needed.\n',kz,VTidx)]);
            end
        end
    end
    if(exist([pwd '/chromthresh/' sprintf('LL5Crvt%d_3afc.mat',VTidx)])==0)
        fprintf('VT measurement for Cr LL5 at JND=%d is still needed.\n',VTidx);
    elseif(sum(bandsign(:))==15)
        findCrVT(VTidx);
        extrafit_Cr(VTidx);
    end
end