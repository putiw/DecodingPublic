clear all; close all; clc

% Dependencies
restoredefaultpath
addpath(genpath([pwd '/helper_functions']));
addpath(genpath('~/Documents/GitHub/TAFKAP')); % https://github.com/Rokers/TAFKAP

bidsDir = '/Volumes/Vision/MRI/DecodingPublic/';

subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID

session = {'01','02','03','04'};
run = [1:10]';
roi =  {'V1','hMT'};
nFolds = 32; % Set to multiple of number of processing cores
maxvoxel = 200;

for iSub = 1:9
    
    sub = subject{iSub};
    
    for iSes = 1:2
        
        ses =  session(iSes*2-1:iSes*2);
        %roi = {'V1','V2','V3','V3A','V3B','hV4','LO1','LO2','hMT','MST','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','VO1','VO2','SPL1','PHC1','PHC2','FEF'};
        params = SetupTAFKAP();
        [allVoxelSample,stim_label] = load_vol(BASE,sub,ses,run,roi,params);
        %
        samples = cell(nFolds,numel(roi,1));
        voxelsize = numel(roi,1);
        
        for i = 1:numel(roi)
            voxelsize(i,1)=size(allVoxelSample{i},2);
            for iF = 1:nFolds
                if voxelsize(i,1) <= maxvoxel
                    samples{iF,i} = allVoxelSample{i};
                else
                    samples{iF,i} = allVoxelSample{i}(:,randperm(voxelsize(i,1),maxvoxel));
                end
            end
        end
        
        if ~exist([bidsDir,'derivatives/dataMat'],'dir')
            mkdir([bidsDir,'derivatives/dataMat'])
        end
        savedata = fullfile(bidsDir,'derivatives/dataMat',[sub,'-ses-' session{iSes*2-1:iSes*2} '-volumn.mat']);
        save(savedata,'samples','stim_label','roi','voxelsize');
        
        %% Run TAFKAP
        
        % Setup design (parameters)
        nScans = length(ses)*length(run); % scans per subject
        nDirs = 8; % motion directions
        params.stimval = stim_label;
        params.runNs = reshape(repmat(1:nScans,nDirs,1),1,[])'; %trainruns; % stimulus block/run
        
        pre = cell(nFolds,1);  % Preallocate
        for ii = 1:nFolds
            p{ii} = params;
            c = cvpartition(params.stimval, 'Holdout', 0.1); % stratify by motion direction, but not scan
            p{ii}.train_trials = c.training;
            p{ii}.test_trials = c.test;
            pre{ii} = params.stimval(c.test);
        end
        
        % Run TAFKAP
        ests = cell(numel(roi),1); % Preallocate
        uncs = cell(numel(roi),1); % Preallocate
        pres = cell(numel(roi),1); % Preallocate
        for whichRoi = 1:numel(roi)
            parfor ii = 1:nFolds
                ii
                rng(ii);% To counter the effects of TAFKAP_Decode setting the system rand seed to const. This was an EXTRAORDINARLY hard bug to find.
                [est{ii}, unc{ii}, liks{ii}, hypers{ii}] = TAFKAP_Decode(samples{ii,whichRoi}, p{ii});
            end
            ests{whichRoi} = cell2mat(est');
            uncs{whichRoi} = cell2mat(unc');
            pres{whichRoi} = cell2mat(pre);
        end
        toc
        %% Calculate performance
        params.subjects = sub;
        saveresult = cell(numel(roi),1);
        for whichRoi = 1:numel(roi)
            % conmat = confusionmat(categorical(ceil(pres{rr}/22.5)),categorical(ceil(ests{rr}/22.5))); % continuous data
            conmat = confusionmat(categorical(pres{whichRoi}),categorical(ests{whichRoi})); % categorical data
            conmat = 100.*conmat./sum(conmat,2); % convert to % accuracy
            conmat = [conmat; conmat(1,:)]; % wrap matrix
            conmat = [conmat, conmat(:,1)];
            saveresult{whichRoi} = conmat';
            disp(['Classification performance ' roi{whichRoi} ': ' num2str(100.*mean(pres{whichRoi}==ests{whichRoi})) '%'])
        end
        if ~exist([bidsDir,'derivatives/resultsMat'],'dir')
            mkdir([bidsDir,'derivatives/resultsMat'])
        end
        f = fullfile(bidsDir,'derivatives/resultsMat',[sub,'-ses-' ses{:} '-TAFKAP.mat']);
        save(f,'saveresult','ests','uncs','pres','roi','voxelsize');
    end
end