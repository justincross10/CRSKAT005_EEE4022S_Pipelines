%% Load WAV file
wavFilePath = 'C:\Users\justm\OneDrive - University of Cape Town\2025\EEE4022S\Testing\Testing 26-09\Gball 4 fast\adc_data.wav';

[raw_data, fs] = audioread(wavFilePath);  % Read the audio data and its sampling frequency from the WAV file
 
raw_data = raw_data(:,1);  % Use only the first channel for processing

%% STFT / Spectrogram Parameters
window = 1024;                              % Define the analysis window length in samples
noverlap = round(0.95*window);              % Define the overlap between consecutive windows
step = window - noverlap;                   % Compute the step size between successive windows
nfft = window;                              % Set the FFT length
w_type = hann(window,'periodic');           % Create a Hann window to reduce spectral leakage in the FFT
num_frames = floor((length(raw_data) - window)/step) + 1; % Calculate how many time frames will fit in the input signal
S = zeros(nfft/2, num_frames);              % Preallocate the spectrogram matrix

%% Compute STFT manually
for k = 1:num_frames
    idx = (k-1)*step + 1;                          % Define the start index of the current frame
    seg = raw_data(idx:idx+window-1) .* w_type;    % Extract the current signal segment and apply the window function
    X = fft(seg, nfft);                            % Compute the FFT of the windowed segment
    P = (abs(X(1:nfft/2)).^2) / sum(w_type.^2);    % Take the magnitude squared of only the first half of the FFT (signal is real valued) and normalize by window energy
    S(:,k) = P;                                    % Store the power spectrum for this frame
end

%% Crop frequency axis
f_lim = 4000;                               % Maximum frequency (Hz) to retain for analysis
freq_axis = (0:(nfft/2-1)) * (fs/nfft);     % Generate the frequency axis for the one-sided FFT
freq_idx = freq_axis <= f_lim;              % Identify indices corresponding to frequencies below the limit
freq_axis_cropped = freq_axis(freq_idx);
P_cropped = 10*log10(S(freq_idx,:) + eps);  % Convert the spectrogram to dB scale for visualization
P_lin_cropped = 10.^(P_cropped / 10);       % Convert back to linear scale if needed for further processing

%% Time axis
time_axis = ((0:num_frames-1)*step + window/2)/fs;

% Time crop
t_start = 0; % start time (s)
t_end = 3;   % end time (s)
time_idx = find(time_axis >= t_start & time_axis <= t_end);

% Apply cropping
time_axis = time_axis(time_idx);
P_cropped = P_cropped(:, time_idx);
P_lin_cropped = P_lin_cropped(:, time_idx);

%% Apply global threshold
alpha = 10;           % Scaling factor for threshold
quiet_frames = 1:5; % Frames with no streak to estimate noise
min_bin = 5;        % Frequency bin cutoff 

% Use only bins above cutoff
bins_to_use = min_bin:size(P_lin_cropped, 1);
noise_est = mean(P_lin_cropped(bins_to_use, quiet_frames), 'all'); % Estimate the average background noise power in the selected frequency bins

% Compute the detection threshold
threshold = alpha * noise_est;
detections = P_lin_cropped > threshold; % Generate a binary detection mask

% Zero-out clutter region
detections(1:5, :) = 0;

%% Plot raw spectrogram
figure;
imagesc(time_axis, freq_axis_cropped, P_cropped);
axis xy; 
colormap parula;
c = colorbar;
c.Label.String = 'Magnitude [dB]';
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Raw Spectrogram of Data');
xlim([t_start t_end]);
dcm = datacursormode(gcf);   % Get data cursor object for current figure
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');  % Show coordinates in a small window

%% Plot the detections after the global threshold has been applied
figure;
imagesc(time_axis, freq_axis_cropped, P_cropped);
axis xy; 
colormap(parula); 
c = colorbar;
c.Label.String = 'Magnitude [dB]';
xlabel('Time (s)'); ylabel('Frequency (Hz)');
xlim([t_start t_end])
title('Spectrogram with Detections from Global Threshold');
hold on;
[f_detect, t_detect] = find(detections);
plot(time_axis(t_detect), freq_axis_cropped(f_detect), 'r.', 'MarkerSize', 8);
hold off;

%% Find connected clusters in the detection map

[rows, cols] = size(detections); % Get the size of the binary detection matrix
visited = zeros(rows, cols);     % Initialize a matrix to keep track of which pixels have been visited
clusters = {};                   % Initialize a cell array to store the clusters found

% Loop over each pixel in the detection map
for r = 1:rows
    for c = 1:cols
        % If the pixel is part of a detection and has not been visited
        if detections(r,c) == 1 && visited(r,c) == 0
            % Start a new cluster using flood-fill
            stack = [r c]; % Stack for depth-first search
            cluster = [];  % Temporary list to store pixels in this cluster

            % Perform flood-fill to find all connected pixels
            while ~isempty(stack)
                [cr, cc] = deal(stack(end,1), stack(end,2));
                stack(end,:) = [];
                if visited(cr, cc) == 0 && detections(cr, cc) == 1          % If the pixel is unvisited and a detection
                    visited(cr, cc) = 1;                                    % Mark it as visited
                    cluster = [cluster; cr, cc];                            % Add it to the current cluster
                    neighbors = [cr-1 cc; cr+1 cc; cr cc-1; cr cc+1];       % Find 4-connected neighbors (up, down, left, right)
                    % Keep only neighbors that are within the matrix bounds
                    neighbors = neighbors(neighbors(:,1)>0 & neighbors(:,1)<=rows & neighbors(:,2)>0 & neighbors(:,2)<=cols, :); 
                    stack = [stack; neighbors];  % Add neighbors to the stack for further exploration
                end
            end
            clusters{end+1} = cluster; % Store the completed cluster in the cell array
        end
    end
end

%% Filter clusters based on minimum size to isolate potential ball targets

minClusterSize = 100;                          % Define the minimum number of pixels required for a cluster to be considered valid
numPixels = cellfun(@(c) size(c,1), clusters); % Count the number of pixels in each detected cluster
validIdx = find(numPixels >= minClusterSize);  % Identify clusters that meet or exceed the minimum size threshold
validClusters = clusters(validIdx);            % Select only the valid clusters for further processing

if isempty(validClusters)
    error('No clusters with at least %d cells found.', minClusterSize); % Give an error if no clusters meet the size requirement
end

%% Plot valid clusters only
figure;
imagesc(time_axis, freq_axis_cropped, detections);
axis xy; 
colormap(flipud(gray));
xlabel('Time (s)'); 
ylabel('Frequency (Hz)');
xlim([t_start t_end]);
title('Detection Map with Valid Clusters for Global Threshold');
hold on;
colors = lines(length(validClusters));

for k = 1:length(validClusters)
    cluster = validClusters{k};
    y = freq_axis_cropped(cluster(:,1));
    x = time_axis(cluster(:,2));
    plot(x, y, '.', 'MarkerSize', 10, 'Color', colors(k,:));
end
hold off;

%% Radar parameters
c = 3e8;             % Speed of light
f_radar = 24.125e9;  % Radar carrier frequency
lambda = c/f_radar;  % Carrier Wavelength
numClusters = length(validClusters);

%% Identify the cluster corresponding to the maximum velocity

allVelocities = zeros(1, numClusters); % Preallocate array to store the peak velocity of each cluster

% Loop through all valid clusters
for k = 1:numClusters
    clusterFreqs = freq_axis_cropped(validClusters{k}(:,1)); % Extract the frequency bins for the current cluster
    [maxFreq, ~] = max(clusterFreqs);                        % Find the maximum frequency within this cluster
    allVelocities(k) = (maxFreq * lambda) / 2;               % Convert the maximum frequency to velocity using the Doppler formula
end

[peak_vel, maxClusterIdx] = max(allVelocities); % Identify the cluster with the highest velocity
predictedVelocity = peak_vel;
fprintf('The predicted ball velocity is %.2f m/s\n', predictedVelocity);

