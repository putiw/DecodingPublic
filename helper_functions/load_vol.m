function [dataset stim_label] = load_vol(BASE,sub,ses,run,roi,p)
%
% Load and process (detrend, normalize) data for classification

% Allocate data
DATA = cell(numel(ses).*numel(run),numel(roi));
stim_label=[];
% true trial design matrix [even run, odd run]
stim = [[5:-1:1 8:-1:6]' [4:8 1:3]'];

for whichSession = 1:numel(ses)
    
    for whichRun = 1:numel(run)
        
        datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{whichSession},'/func/', ...
            sub,'_ses-',ses{whichSession},'_task-3dmotion_run-',num2str(whichRun), ...
            '_space-T1w_desc-preproc_bold.nii.gz'];
        
        disp(['Loading: ' datapath]);
        Func = niftiread(fullfile(datapath));
        
        
        % Drop initial frames to eliminate transients and reach steady state
        framesToDrop = 10;
        Func = Func(:,:,:,framesToDrop+1:end); % Drop n frames
        numFrames = size(Func,4);
        
        % Extract timecourses within the ROIs
        for whichRoi = 1:numel(roi)
            
            roiPath = [BASE,'derivatives/fmriprep/',sub,'/ses-01/anat/rois/', ...
                sub,'_space-T1w_downsampled_',roi{whichRoi},'.nii.gz'];
            
            ROI = niftiread(fullfile(roiPath));
            roiSize = length(find(ROI));
            [x y z] = ind2sub(size(ROI),find(ROI));
            
            % Extract raw intensity timeseries
            samples = zeros(numFrames,roiSize);
            for voxel = 1:roiSize
                samples(:,voxel) = squeeze(Func(x(voxel),y(voxel),z(voxel),:));
            end
            
            %%  detrend + normalize
            %samples = detrend(samples,1); % linear detrend
            
            % fft-based detrend (works better than linear detrend)
            fmriFFT = fft(samples);
            fmriFFT(1:5,:) = zeros(5,roiSize);
            fmriFFT(end-4:end,:) = zeros(5,roiSize);
            samples = real(ifft(fmriFFT));
            
            % samples = samples-mean(samples,2); % ROI-average based detrend
            % Used in TAFKAP. Does not do much beyond linear detrend for
            % classify
            
            samples = (samples(1:2:end-1,:) + samples(2:2:end,:)) ./2; % take average of every 2 TRs
            
            samples = normalize(samples); % z-score samples
            %                     % needed for TAFKAP, does not do much for classify
            %
            % Organize as scan-based or trial-based analysis
            switch p.sample_unit
                case {'scan'}
                    samples = squeeze(mean(reshape(samples,8,15,[]),2)); % average every 8th datapoint
                case {'trial'}
                    % do nothing, keep samples
                otherwise
                    error('Unknown sample unit')
            end
            
            %  samples = normalize(samples); % z-score samples
            % needed for TAFKAP, does not do much for classify
            
            % output is downsampled from TRs (240) x voxels to average response per run per direction (8) x voxels
            
            DATA{numel(run)*(whichSession-1)+whichRun,whichRoi} = samples;
        end % end of roi
        % end of run
        
        % Return stimulus labels
        stim_label =  [stim_label; stim(:,mod(whichRun,2)+1)]; % assign stimulus labels based on odd/even run        
        
    end
    
end

dataset = cell(1,numel(roi));
for whichRoi = 1:numel(roi)
    dataset{whichRoi} = cell2mat(DATA(:,whichRoi));
end

end