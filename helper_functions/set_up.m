
function set_up(bidsDir,gitDir,freesurferDir)

% github

addpath(genpath([gitDir '/gifti']))
addpath(genpath(fullfile(gitDir, 'knkutils'))); % https://github.com/cvnlab/knkutils
addpath(genpath(fullfile(gitDir, 'cvncode'))); % https://github.com/cvnlab/cvncode

% freesurfer settings
fslDir = '/usr/local/fsl';
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fslDir '/bin']); % add freesurfer/bin to path
setenv('PATH', [PATH ':' fslDir '/bin']); % add freesurfer/bin to path
setenv('PATH', sprintf('/usr/local/bin:%s', getenv('PATH'))); % add /usr/local/bin to PATH
setenv('PATH', ['/Applications/freesurfer/7.2.0/bin:' getenv('PATH')]);
setenv('PATH', [PATH ':/usr/local/fsl/bin']);
setenv('FSLDIR', fslDir);
setenv('FREESURFER_HOME', freesurferDir);
addpath(genpath(fullfile(freesurferDir, 'matlab')));
setenv('SUBJECTS_DIR', [bidsDir '/derivatives/freesurfer']);

end