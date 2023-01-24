% Extract probabilistic ROIs
% Before running this script, ensure that fmriprep has been run
% Requires matlab is started from the terminal
% /Applications/MATLAB_R2017b.app/bin/matlab or equivalent
%
% Dependencies: Freesurfer, FSL
%
% Setup
setenv('PATH','$PATH:/bin:/usr/bin:/Applications/freesurfer/bin:/usr/local/fsl/bin:/Applications/freesurfer/7.2.0/bin/');
setenv('FSL_DIR', '/usr/local/fsl');
setenv('FREESURFER_HOME','/Applications/freesurfer/7.2.0');  % this to tell where FSL folder is
clear all

% Terminal commands for freesurfer install
% export FREESURFER_HOME=/Applications/freesurfer/7.2.0
% export SUBJECTS_DIR=$FREESURFER_HOME/subjects
% source $FREESURFER_HOME/SetUpFreeSurfer.sh

%% User-changable variables
projectDir      = '/Volumes/Vision/MRI/Decoding';
subjectName = '0201';
sessionName = '01';
spaceName   = 'T1w';

%% Make no changes below here
sub = 'sub-0201'; %subject's freesurfer directory name
ses = ['sub-' subjectName];
space = ['space-' spaceName];
anatDir = fullfile(projectDir, '/derivatives/fmriprep/', sub, ses, '/anat'); %anatomy directory with subfolder /ROIs for shared ROIs
roiDir = [anatDir '/rois'];  % ADD A CHECK TO SEE IF THIS EXISTS, AND IF NOT, CREATE IT.

%% Wang atlas
%
% Based on https://scholar.princeton.edu/napl/resources

roi = {'V1v';'V1d';'V2v';'V2d';'V3v';'V3d';'hV4';'VO1';'VO2';'PHC1';'PHC2'; ...
    'MST';'hMT';'LO2';'LO1';'V3B';'V3A';'IPS0';'IPS1';'IPS2';'IPS3';'IPS4'; ...
    'IPS5';'SPL1';'FEF'};
index = [1:length(roi)]';
wangLUT = table(index, roi);

%% Save atlas as individual rois

system(['/usr/local/bin/docker run -ti --rm -v ' ...
    projectDir '/derivatives/freesurfer:/subjects' ...
    ' nben/neuropythy' ...
    ' atlas --verbose ' sub ' --volume-export']);

%% Create rois and match resolution in Project
%
% Instead of downsampled consider adding ref-anat and ref-func to the
% filenames
if ~exist(roiDir, 'dir')
    mkdir(roiDir);
end

for ii = 1:length(wangLUT.index)
    
    fsRoiDir = fullfile(projectDir, 'derivatives/freesurfer', sub, 'mri/rois');
    inFile = fullfile(projectDir, 'derivatives/freesurfer', sub, 'mri/wang15_mplbl.mgz');
    roiFile = fullfile(roiDir, [sub '_' space '_' wangLUT.roi{ii} '.nii.gz']); %fullfile(fsRoiDir,[wangLUT.roi{ii} '.mgz']);
    refFile = fullfile(anatDir, [sub '_' ses '_acq-highres' '_desc-preproc_T1w.nii.gz']); % use space ref
    system(['mri_binarize --i ' inFile ...
        ' --o ' roiFile ...
        ' --match ' num2str(wangLUT.index(ii))]);
    
    % Doing a straight up mri_convert on T1w space with the func as ref segfaults, so
    % workaround
    
    outFile = fullfile(roiDir, [sub '_' space '_' wangLUT.roi{ii} '.nii.gz']);
    system(['mri_convert -rl ' refFile ' -rt nearest ' roiFile ' ' outFile]);
    
    % Downsample to functional resolution - TODO: Change downsampled to ref-anat
    % and ref-func
    refFile = dir(fullfile(projectDir, 'derivatives/fmriprep', sub, ses, 'func', [sub '_' ses '_task-*_run-1_' space '_boldref.nii.gz'])); % uses fmriprep coreg ref
    refFile = fullfile(refFile.folder, refFile.name);
    roiFile = outFile;
    outFile = fullfile(roiDir, [sub '_' space '_downsampled_' wangLUT.roi{ii} '.nii.gz']);
    system(['mri_convert -rl ' refFile ' -rt nearest ' roiFile ' ' outFile]);
end

% Combine ROIs
roiC = {'V1', 'V2', 'V3'};
for ii = 1:3
    roi1 = fullfile(roiDir, [sub '_' space '_downsampled_' wangLUT.roi{2*ii-1} '.nii.gz']);
    roi2 = fullfile(roiDir, [sub '_' space '_downsampled_' wangLUT.roi{2*ii} '.nii.gz']);
    outFile = fullfile(roiDir, [sub '_' space '_downsampled_' roiC{ii} '.nii.gz']);
    
    system(['fslmaths ' roi1 ' -add ' roi2 ' ' outFile]);
end
