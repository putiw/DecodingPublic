%% set up
clear all; close all; clc
restoredefaultpath;
addpath(genpath([pwd '/helper_functions']));
addpath(genpath('~/Documents/GitHub/TAFKAP')); % https://github.com/Rokers/TAFKAP
bidsDir = '/Volumes/Vision/MRI/DecodingPublic/';
subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};    % subject ID
session = {'01','02','03','04'};
run = [1:10]';
roi = {'V1','V2','V3','V3A','V3B','hV4','LO1','LO2','hMT','MST','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','VO1','VO2','SPL1','PHC1','PHC2','FEF'};
nFolds = 32; % Set to multiple of number of processing cores
maxvoxel = inf; %max voxel used in the analysis
if ~exist([bidsDir,'derivatives/dataMat'],'dir')
    mkdir([bidsDir,'derivatives/dataMat'])
end
if ~exist([bidsDir,'derivatives/resultsMat'],'dir')
    mkdir([bidsDir,'derivatives/resultsMat'])
end
%% TAFKAP
for iSub = 1:9
    for iSes = 1:2
        ses =  session(iSes*2-1:iSes*2); % session ID        
        
        % load data
        dataFile = fullfile(bidsDir,'derivatives/dataMat',[subject{iSub},'-ses-' session{iSes*2-1:iSes*2} '-volumn.mat']);
        if ~exist(dataFile,'file') % if cleaned .mat file doesn't exist, load from gifti file
            [allVoxelSample,stim_label] = load_vol(bidsDir,subject{iSub},ses,run,roi);
            disp(['saving ' subject{iSub},'-ses-' session{iSes*2-1:iSes*2} '-volumn.mat to ' bidsDir,'/derivatives/dataMat/'])
            save(dataFile,'allVoxelSample','stim_label','roi');
        else  % load cleaned .mat file if it exists to save time
            disp(['data file exists, loading ' subject{iSub},'-ses-' session{iSes*2-1:iSes*2} '-volumn.mat from ' bidsDir,'/derivatives/dataMat/'])
            load(dataFile);
        end
        
        % load data if using diffferent num of voxels
        samples = cell(nFolds,numel(roi,1));
        voxelsize = zeros(numel(roi),1);
        % random select *maxvoxel* number voxels from each ROIs for each of the
        % cross-validatin fold
        for iRoi = 1:numel(roi)
            voxelsize(iRoi,1)=size(allVoxelSample{iRoi},2);
            for iF = 1:nFolds
                if voxelsize(iRoi,1) <= maxvoxel % if a given ROI has less than *maxvoxel* number of voxels, then use all the entire ROI
                    samples{iF,iRoi} = allVoxelSample{iRoi};
                else
                    samples{iF,iRoi} = allVoxelSample{iRoi}(:,randperm(voxelsize(iRoi,1),maxvoxel));
                end
            end
        end
        
        %% Run TAFKAP
        
        resultFile = fullfile(bidsDir,'derivatives/resultsMat',[subject{iSub},'-ses-' ses{:} '-TAFKAP.mat']);
        if exist(resultFile,'file')
            prompt = sprintf('This file already exists:\n%s\nDo you want to overwrite it?', resultFile);
            question = 'Overwrite?';
            answer = questdlg(prompt, question, 'Yes','No','Yes');
            switch answer
                case 'yes' % re-analysis the data and overwrite the results
                    % Setup design (parameters)
                    params = SetupTAFKAP(); % set up parameters
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
                    % Calculate performance
                    params.subjects = subject{iSub};
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
                    
                    save(resultFile,'saveresult','ests','uncs','pres','roi','voxelsize');
                    
                case 'No' % result file already exists, print results
                    load(resultFile);
                    for whichRoi = 1:numel(roi)
                        disp(['Classification performance ' roi{whichRoi} ': ' num2str(100.*mean(pres{whichRoi}==ests{whichRoi})) '%'])
                    end
            end
        end
    end
end