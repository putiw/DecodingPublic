clear all;clc;close all; restoredefaultpath;
addpath(genpath([pwd '/helper_functions']));
bidsDir = '/Volumes/Vision/MRI/DecodingPublic';
fsDir = '/Applications/freesurfer/7.2.0';
gitDir = '~/Documents/GitHub';
set_up(bidsDir,gitDir,fsDir)

tic
%% find searchlight index

% load fsaverage6 .inflated surf
surfL = fullfile([bidsDir '/derivatives/freesurfer/fsaverage6/surf/lh.inflated']);
surfL = fs_read_surf(surfL);
surfR = fullfile([bidsDir '/derivatives/freesurfer/fsaverage6/surf/rh.inflated']);
surfR = fs_read_surf(surfR);
num = size(surfL.coord,2);
% find what is inside for each searchlight
idx = zeros(num,100,2);
for iVer = 1:num
    distance = sqrt(sum((surfL.coord - surfL.coord(:,iVer)).^2));
    [~,idx(iVer,:,1)] = mink(distance,100);
    distance = sqrt(sum((surfR.coord - surfR.coord(:,iVer)).^2));
    [~,idx(iVer,:,2)] = mink(distance,100);
end

%%
subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};     % subject ID
ses =   {'01','02'}; %{'03','04 '}
run = [1:10]';
hemi = {'L','R'};

for iSub = 1:9
    sub = subject{iSub};
    
    %% load surf data
    
    samples = load_surf(bidsDir,sub,ses,run,'fsaverage6',hemi,'fft');
    label = [repmat([[4:8 1:3]';[5:-1:1 8:-1:6]'],10,1) repmat([8 1:7 1 8:-1:2]',10,1)];
    
    if ~exist([bidsDir,'/derivatives/dataMat'],'dir')
        mkdir([bidsDir,'/derivatives/dataMat']);
    end
    savedata = fullfile(bidsDir,'/derivatives/dataMat',[sub,'-ses-' ses{:} '-fsaverage6.mat']);
    save(savedata,'samples','label');
    
    %% MATLAB Classify
    
    nRep = 300; % bootstraps repeats
    test_n = 1; % percetange of data set as testing trials
    
    trial_n = 160;  % number of trials
    dir_n = 8; % number of directions
    nScans = trial_n/dir_n; % number of scans
    
    Training_group = sort(repmat(repelem(1:8,1)',nScans-test_n,1)); % training dataset stimulus label
    Testing_group = repmat(repelem(1:8,1)',test_n,1); % testing dataset stimulus label
    
    pc = zeros(nRep,num,2); % percentage correct
    
    for iHemi = 1:2
        
        parfor ii = 1:num
            
            data = samples(:,idx(ii,:,iHemi),iHemi); % 160 (trials) x 100 (vertices)
            
            for kk = 1:nRep
                
                Training_data = [];
                Testing_data = [];
                
                for mm = 1:dir_n
                    train_run = 1:nScans; % run index
                    test_run = randperm(nScans,test_n); %random select (test_n) trial for each direction as testing trial index (without replacement)
                    train_run(test_run) = []; % set the rest as the training trial index
                    dir_idx = find(label(:,iHemi)==mm); % find trial index for this direction
                    Testing_data = [Testing_data; data(dir_idx(test_run),:)];
                    Training_data = [Training_data; data(dir_idx(train_run),:)];
                end
                
                decoded_group = classify(Testing_data,Training_data,Training_group,'diaglinear');
                pc(kk,ii,iHemi) = length(find(decoded_group == Testing_group))/length(Testing_group); % percentage correct
                [sub iHemi ii]
            end
        end
        
    end
    
    acc = nanmean(pc,1);
    acc = cat(1,acc(:));
    
    if ~exist([bidsDir,'/derivatives/resultMat'],'dir')
        mkdir([bidsDir,'/derivatives/resultMat']);
    end
    savedata = fullfile(bidsDir,'/derivatives/resultMat',[sub,'_ses-' ses{:} '_searchlight_fsaverage6.mat']);
    save(savedata,'acc');
    
    
    
end
toc