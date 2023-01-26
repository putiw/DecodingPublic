function p = SetupTAFKAP()

% Default settings for parameters in 'p'
p.Nboot = 500; % Maximum number of bootstrap iterations
% For discrete number of stimulus categories, two basis functions
% cannot be equidistant from a single category (so choose 1 or 5 in
% our case)
p.precomp_C = 1; %4; % How many sets of channel basis functions to use (swap between at random) - for PRINCE, this value is irrelevant
p.randseed = 1234; % The seed for the (pseudo-)random number generator, which allows the algorithm to reproduce identical results whenever it's run with the same input, despite being stochastic.
p.prev_C = false; % Regress out contribution of previous stimulus to current-trial voxel responses?
p.dec_type = 'TAFKAP'; % 'TAFKAP' or 'PRINCE'
p.DJS_tol = 1e-8; % If the Jensen-Shannon Divergence between the new likelihoods and the previous values is smaller than this number, we stop collecting bootstrap samples (before the maximum of Nboot is reached). If you don't want to allow this early termination, you can set this parameter to a negative value.
p.nchan = 2; % Number of "channels" i.e. orientation basis functions used to fit voxel tuning curves
p.chan_exp = 5; % Exponent to which basis functions are raised (higher = narrower)

p.stim_type = 'categorical'; % values in p.stim_val must be integers corresponding to category labels, data is treated as continuous if omitted

% p.stimval = ds.sa.targets*180/8;
% p.runNs = ds.sa.chunks;

% Scan-based or trial-based decoding
% Trial-based decoding works (using fft-based detrending), but gives
% overall poorer decoding performance (~20%) vs (~60%)
% raw returns data in x,y,z,t format
p.sample_unit = 'trial'; %'trial'; % 'scan'

end