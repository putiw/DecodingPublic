function DATA = load_surf(dataDir,sub,ses,run,scantype,hemi,detrendcase)

% Inputs:
%     dataDir:     Path to project folder containing the derivative folder; E.g. '/Volumes/Vision/MRI/Decoding'
%     sub:            Subject ID; E.g. 'sub-0201'
%     ses:            Session ID; E.g. {'01','02'}
%     run:            Run ID; E.g. [1:10]'
%     scantype:       Data space; E.g. 'fsnative', 'fsaverage6', etc
%     hemi:           Hemisphere; E.g. {'L','R'}
%     detrendcase:    Different detrending options; E.g. 'fft', 'linear', 'roi-average', etc

% Outputs:
%     samples:        data (8 directions by 20 runs by 40962 vertices by 2 hemisphere)

% Allocate data
DATA_temp = cell(numel(ses).*numel(run),numel(hemi));
DATA = zeros(8*numel(ses).*numel(run),40962,numel(hemi));
for whichHemi = 1:numel(hemi)
    for whichSession = 1:numel(ses)
        for whichRun = 1:numel(run)
            disp(['Loading: ' sub '-ses-' ses{whichSession} '-run-' num2str(whichRun) '-' hemi{whichHemi}])
            
            datapath  = fullfile(dataDir, ['derivatives/fmriprep/', sub, '/ses-', ses{whichSession}, '/func/', ...
                sub,'_ses-', ses{whichSession}, '_task-3dmotion_run-', num2str(whichRun), ...
                '_space-', scantype,'_hemi-', hemi{whichHemi} '_bold.func.gii']);
            
            Func  = gifti(datapath); Func  = Func.cdata;
            
            % Drop initial frames to eliminate transients and reach steady state
            framesToDrop = 10;
            samples  = Func(:,framesToDrop+1:end)'; % Drop n frames
            [numFrame, roiSize] = size(samples);
            
            %%  detrend + normalize
            
            switch(detrendcase)
                
                case 'fft'
                    % fft-based detrend (works better than linear detrend)
                    fmriFFT = fft(samples);
                    fmriFFT(1:5,:) = zeros(5,roiSize);
                    fmriFFT(end-4:end,:) = zeros(5,roiSize);
                    samples = real(ifft(fmriFFT));
                    
                case 'linear'
                    % linear detrend
                    samples = detrend(samples,1);
                    
                case 'roi-average'
                    % ROI-average based detrend used in TAFKAP.
                    % Does not do much beyond linear detrend for classify
                    samples = samples-mean(samples,2);
                    
            end
            
            samples = (samples(1:2:end-1,:) + samples(2:2:end,:)) ./2; % take average of every 2 TRs
            samples = normalize(samples); % z-score sample
            % needed for TAFKAP, does not do much for classify
            
            samples = squeeze(mean(reshape(samples,8,15,[]),2)); % average every 8th datapoint
            % output is downsampled from TRs (240) x voxels to average response per run per direction (8) x voxels
            DATA_temp{numel(run)*(whichSession-1)+whichRun,whichHemi} = samples;
            
        end
    end
    DATA(:,:,whichHemi)  = cell2mat(DATA_temp(:,whichHemi));
    
end
end