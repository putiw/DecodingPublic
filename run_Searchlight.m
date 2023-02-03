%% set up

clear all;clc;close all;
restoredefaultpath;
addpath(genpath([pwd '/helper_functions'])); % path to helper functions

% location variables cannot have ~ in them
bidsDir = '/Volumes/Vision/MRI/DecodingPublic'; % '/Users/rokers/Documents/OpenNeuro/ds004443-download'; % path to bids
fsDir = '/Applications/freesurfer/7.2.0'; % path to freesurfer
gitDir = '/Users/rokers/Documents/GitHub'; % path to github
set_up(bidsDir,gitDir,fsDir) % set up path and dependencies

%% find searchlight index

% load fsaverage6 .inflated surf to get xyz coordinates
surfL = fullfile([bidsDir '/derivatives/freesurfer/fsaverage6/surf/lh.inflated']);
surfL = fs_read_surf(surfL);
surfR = fullfile([bidsDir '/derivatives/freesurfer/fsaverage6/surf/rh.inflated']);
surfR = fs_read_surf(surfR);
num = size(surfL.coord,2);

% find what is inside for each searchlight using the xyz coordinates
idx = zeros(num,100,2); %for each vertex, list the closest 100 vertices, for each hemi
for iVer = 1:num
    distance = sqrt(sum((surfL.coord - surfL.coord(:,iVer)).^2));
    [~,idx(iVer,:,1)] = mink(distance,100);
    distance = sqrt(sum((surfR.coord - surfR.coord(:,iVer)).^2));
    [~,idx(iVer,:,2)] = mink(distance,100);
end

%% Run searchlight

% mkdir
if ~exist([bidsDir,'/derivatives/dataMat'],'dir')
    mkdir([bidsDir,'/derivatives/dataMat']);
end
if ~exist([bidsDir,'/derivatives/resultMat'],'dir')
    mkdir([bidsDir,'/derivatives/resultMat']);
end

% set up params
subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};     % subject ID
session = {'01','02','03','04'}; % which session
run = [1:10]'; % which run
hemi = {'L','R'}; % which hemi
nRep = 1000; % bootstrap repeats
nTest = 1; % percentage of data set as testing trials
nTrial = 160;  % number of trials
nDir = 8; % number of directions
nScans = nTrial/nDir; % number of scans
labelTrain = sort(repmat(repelem(1:8,1)',nScans-nTest,1)); % training dataset stimulus label
labelTest = repmat(repelem(1:8,1)',nTest,1); % testing dataset stimulus label


for iSub = 1:numel(subject) % looping through subjects

    for iSes = 1:2 % horizontal conditions, then vertical conditions

        ses =  session(iSes*2-1:iSes*2);

        % load surf data
        dataFile = fullfile(bidsDir,'/derivatives/dataMat',[subject{iSub},'-ses-' ses{:} '-fsaverage6.mat']);

        if ~exist(dataFile,'file') % if cleaned .mat file doesn't exist, load from gifti file

            % samples -> 160 (trials) x 40962 (vertices) x 2 (hemi)
            samples = load_surf(bidsDir,subject{iSub},ses,run,'fsaverage6',hemi,'fft');
            % stimulus label -> 160 (trials) x 2 (hemi)
            label = [repmat([[4:8 1:3]';[5:-1:1 8:-1:6]'],10,1) repmat([8 1:7 1 8:-1:2]',10,1)];

            disp(['saving ' subject{iSub},'-ses-' ses{:} '-fsaverage6.mat to ' bidsDir,'/derivatives/dataMat/'])
            save(dataFile,'samples','label');

        else  % load cleaned .mat file if it exists to save time
            disp(['data file exists, loading ' subject{iSub},'-ses-' ses{:} '-fsaverage6.mat from ' bidsDir,'/derivatives/dataMat/'])
            load(dataFile);
        end


        % MATLAB Classify
        resultFile = fullfile(bidsDir,'/derivatives/dataMat',[subject{iSub},'-ses-' ses{:} '-fsaverage6.mat']);

        if exist(resultFile,'file')
            prompt = sprintf('This file already exists:\n%s\nDo you want to overwrite it?', resultFile);
            question = 'Overwrite?';
            answer = questdlg(prompt, question, 'Yes','No','Stop','Yes');
            switch answer
                case 'Yes' % re-analysis the data and overwrite the results
                    pc = zeros(nRep,num,2); % percentage correct. bootsrap repeats x vertices x 2 hemi
                    for iHemi = 1:2
                        parfor ii = 1:num

                            data = samples(:,idx(ii,:,iHemi),iHemi); % 160 (trials) x 100 (vertices)

                            for kk = 1:nRep

                                dataTrain = [];
                                dataTest = [];

                                for mm = 1:nDir
                                    train_run = 1:nScans; % run index
                                    test_run = randperm(nScans,nTest); %random select (nTest) trial for each direction as testing trial index (without replacement)
                                    train_run(test_run) = []; % set the rest as the training trial index
                                    whichTrial = find(label(:,iHemi)==mm); % find trial index for this direction
                                    dataTest = [dataTest; data(whichTrial(test_run),:)];
                                    dataTrain = [dataTrain; data(whichTrial(train_run),:)];
                                end

                                decoded_group = classify(dataTest,dataTrain,labelTrain,'diaglinear');
                                pc(kk,ii,iHemi) = length(find(decoded_group == labelTest))/length(labelTest); % percentage correct
                            end
                        end
                    end

                    acc = nanmean(pc,1);
                    acc = cat(1,acc(:));
                    save(resultFile,'acc');

                case 'No' % result file already exists, skip this subject
                case 'Stop'
                    return
            end
        end
    end
end
