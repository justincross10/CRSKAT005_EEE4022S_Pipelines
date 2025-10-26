%% Load 4-channel WAV file and extract specific channels
wavFilePath = 'C:\Users\justm\OneDrive - University of Cape Town\2025\EEE4022S\Testing\Testing 03-10\Radar on bin\Rball 4 1.3m\adc_data.wav';

[raw_data, fs] = audioread(wavFilePath); % Read the audio data and sampling frequency

% Extract individual radar channels for processing
I1 = raw_data(:,4);  
I2 = raw_data(:,1);

%% STFT / Spectrogram Parameters
window = 2048;                              % Define the analysis window length in samples
noverlap = round(0.95*window);              % Define the overlap between consecutive windows
step = window - noverlap;                   % Compute the step size (hop length) between successive windows
nfft = window;                              % Set the FFT length
w_type = hann(window,'periodic');           % Create a Hann window to reduce spectral leakage in the FFT
num_frames = floor((length(I1) - window)/step) + 1;  % Calculate how many time frames will fit in the input signal

S1 = zeros(nfft, num_frames);               % Preallocate the spectrogram matrix for channel I1
S2 = zeros(nfft, num_frames);               % Preallocate the spectrogram matrix for channel I2

%% Compute STFT for I-only signals
for k = 1:num_frames
    idx = (k-1)*step + 1; % Compute the starting index of the current time frame

    % Extract the current windowed segment for each I-channel
    seg1 = I1(idx:idx+window-1) .* w_type;
    seg2 = I2(idx:idx+window-1) .* w_type;

    % Compute FFT of each windowed segment
    X1 = fft(seg1, nfft);
    X2 = fft(seg2, nfft);

    % Shift zero frequency component to the center of the spectrum
    % and store in spectrogram matrices
    S1(:,k) = fftshift(X1);  
    S2(:,k) = fftshift(X2);
end

%% Crop frequency axis to region of interest

f_min = 0;        % Minimum frequency to retain (Hz)
f_max = 5000;     % Maximum frequency to retain (Hz)

freq_axis = (-nfft/2:nfft/2-1)*(fs/nfft);               % Generate the full frequency axis for the FFT after fftshift
freq_idx = (freq_axis >= f_min) & (freq_axis <= f_max); % Identify indices corresponding to the desired frequency range
freq_axis_cropped = freq_axis(freq_idx);                % Keep only the frequencies within the specified range

% Crop the spectrograms of each I-channel to the frequency range
S1 = S1(freq_idx,:);
S2 = S2(freq_idx,:);

P_cropped = 20*log10(abs(S1) + eps);  % Convert the spectrogram to dB scale for visualization
P_lin_cropped = abs(S1).^2; % Compute the linear power spectrogram from the FFT amplitudes
% Only I1 channel used for velocity information

%% Crop time axis to region of interest

time_axis = ((0:num_frames-1)*step + window/2)/fs; % Generate the time axis corresponding to the center of each STFT window
time_idx = (time_axis >= 0 & time_axis <= 3);      % Select only the time frames within the desired range
time_axis = time_axis(time_idx);                   % Keep only the selected time frames in the time axis

% Crop the spectrograms of both I-channels to the selected time range
S1 = S1(:,time_idx);
S2 = S2(:,time_idx);

P_cropped = P_cropped(:, time_idx); 
P_lin_cropped = P_lin_cropped(:,time_idx);

%% 1D Cell-Averaging CFAR along time axis
% CFAR parameters
num_guard = 400;      % Number of guard frames on each side of the Cell Under Test (CUT)
num_ref   = 125;      % Number of reference frames on each side used to estimate noise
alpha     = 50;      % Scaling factor for noise estimate

[num_bins, num_frames] = size(P_lin_cropped); % Get size of the spectrogram
detections = zeros(num_bins, num_frames);     % Preallocate binary detection matrix (1 indicates detection, 0 indicates no detection)

% Loop over frequency bins
for f = 1:num_bins
    for t = (num_ref+num_guard+1):(num_frames-num_ref-num_guard) % Loop over time frames, avoiding edges where reference/guard cells don't fit
        t_range = (t-num_ref-num_guard):(t+num_ref+num_guard);   % Define the full window around the CUT
        
        % Remove guard cells and the CUT itself from reference cells
        ref_frames = t_range;
        ref_frames(num_ref+1:num_ref+2*num_guard+1) = []; 
        
        % Compute the local noise estimate using the reference cells only)
        noise_local = mean(P_lin_cropped(f, ref_frames));
        
        if P_lin_cropped(f,t) > alpha * noise_local  % Apply CFAR detection (mark CUT as a detection if it exceeds threshold)
            detections(f,t) = 1;
        end
    end
end

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
dcm = datacursormode(gcf);   % Get data cursor object for current figure
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');  % Show coordinates in a small window

%% Plot CFAR detections overlay
figure;
imagesc(time_axis, freq_axis_cropped, P_cropped);
axis xy; 
colormap(parula); 
colorbar;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)');
title('Spectrogram with CFAR Detections');
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
minClusterSize = 300;                          % Define the minimum number of pixels required for a cluster to be considered valid
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
title('Detection Map with Valid Clusters');
hold on;
colors = lines(length(validClusters));

for k = 1:length(validClusters)
    cluster = validClusters{k};
    y = freq_axis_cropped(cluster(:,1));
    x = time_axis(cluster(:,2));
    plot(x, y, '.', 'MarkerSize', 10, 'Color', colors(k,:));
end
hold off;

%% Identify the cluster with the highest frequency detection
maxClusterFreqs = cellfun(@(c) max(freq_axis_cropped(c(:,1))), validClusters); % Max frequency in each cluster
[~, highestFreqClusterIdx] = max(maxClusterFreqs);                             % Index of cluster with highest frequency
highestFreqCluster = validClusters{highestFreqClusterIdx};                     % Extract that cluster

%% Radar parameters
c = 3e8;             % Speed of light
f_radar = 24.125e9;  % Radar carrier frequency
lambda = c/f_radar;  % Carrier Wavelength
d = 0.008763;        % Distance between KLC7 receivers (8.763mm)

%% Precompute analytic signals for AoA calculation
% Copy the original I-channel signals for Angle-of-Arrival (AoA) processing
I1_aoa = I1;  
I2_aoa = I2;

%% Determine frequency range of the selected cluster
freq_bins = highestFreqCluster(:,1);             % Extract the frequency bin indices
freq_values = freq_axis_cropped(freq_bins);      % Convert to actual frequency values in Hz
maxFreqCluster = max(freq_values);               % Maximum frequency in this cluster

%% Compute maximum velocity based on the maximum frequency bin in the cluster
predictedVelocity = (maxFreqCluster * lambda) / 2;  % Velocity equation using Doppler shifts

%% Design and apply bandpass Butterworth filter
% Define adaptive bandpass limits around the cluster based off where
% maximum frequency bin occurs
f_low = max(maxFreqCluster - 1500, 0);             % Lower threshold
f_high = maxFreqCluster + 100;                    % Upper threshold

fprintf('Adaptive Bandpass: %.1f Hz – %.1f Hz\n', f_low, f_high);

[b, a] = butter(4, [f_low f_high] / (fs/2), 'bandpass'); % Design a 4th-order bandpass Butterworth filter

% Apply zero-phase filtering to both I-channels
I1_aoa = filtfilt(b,a,I1_aoa);
I2_aoa = filtfilt(b,a,I2_aoa);

%% Create analytic signals using Hilbert transform
s1 = hilbert(I1_aoa);
s2 = hilbert(I2_aoa);

%% Compute instantaneous AoA using strongest CFAR pixel per frame

freq_thresh = f_low; % Minimum frequency threshold to consider

unique_times = unique(highestFreqCluster(:,2));  % Extract unique time frames from the hgihest frequency cluster

% Initialize arrays to store results
AoA_all = [];  % Stores [time_sec, freq_Hz, AoA_deg] for each frame

% Loop over each time frame in the cluster
for k = 1:length(unique_times)
    t_bin = unique_times(k);
    pixels = highestFreqCluster(highestFreqCluster(:,2)==t_bin,1);     % Get all pixels in this cluster at this time frame

    % Only keep bins above the frequency threshold
    pixels = pixels(freq_axis_cropped(pixels) >= freq_thresh);
    if isempty(pixels)
        continue; % Skip frame if no pixels meet threshold
    end
    
    % Map STFT frame to corresponding I-sample index (center of window)
    idx_center = (t_bin-1)*step + round(window/2);
    if idx_center < 1 || idx_center > length(I1)
        continue; % Skip if index is out of bounds
    end
    
    % Compute instantaneous phase difference between channels
    C = s2(idx_center) * conj(s1(idx_center));
    dphi = angle(C);
    arg = (dphi*lambda)/(2*pi*d); % Convert phase difference to AoA (radians) using antenna spacing and wavelength

    % Only compute AoA if argument is within valid range [-1, 1]
    if abs(arg) <= 1
        theta_deg = asin(arg)*180/pi; % Convert to degrees
        AoA_all = [AoA_all; time_axis(t_bin), theta_deg]; %, freq_axis_cropped(f_bin)
    end
end

%% Plot instantaneous AoA vs time
figure; 
hold on; 
grid on;
plot(AoA_all(:,1), AoA_all(:,2), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor','b');
xlabel('Time (s)');
ylabel('Elevation AoA [deg]');
title('Instantaneous AoA vs Time');
ylim([-20 20]); % adjust according to expected angle range

%% Smooth AoA over time using a moving average
if ~isempty(AoA_all)
    AoA_smooth = AoA_all(:,2); % Extract the raw AoA measurements from the results

    % Apply a moving average filter to reduce noise and small fluctuations
    smooth_window = 21;                              % Number of frames over which to average
    AoA_smooth = movmean(AoA_smooth, smooth_window);
end

%% Smooth AoA over time using Kalman filter and RTS smoother

% Extract time and measured AoA from the results
t = AoA_all(:,1);        % Time vector (seconds)
AoA_meas = AoA_all(:,2);        % Measured AoA (degrees)

if ~isempty(AoA_all)

    % Estimate noise levels for Kalman filter             
    dAoA = diff(AoA_meas);
    R = var(dAoA);                   % Measurement noise covariance
    Q = R * 0.01;                    % Process noise covariance

    % Robust initialization
    n = length(AoA_meas);
    win = min(5, n);                            % Number of initial samples to compute median
    init_val = median(AoA_meas(1:win));         % Robust initial AoA estimate
    init_var = max(var(AoA_meas(1:win)), 1e-3); % Initial variance to avoid zero

    % Preallocate arrays for filter and smoother
    x_filt = zeros(1, n);   % Filtered state
    P_filt = zeros(1, n);   % Filtered covariance
    x_pred = zeros(1, n);   % Predicted state
    P_pred = zeros(1, n);   % Predicted covariance

    % Initialize first state with robust estimate
    x_filt(1) = init_val;            
    P_filt(1) = init_var * 100; % Large initial uncertainty for flexibility

    % Forward pass (Kalman filter)
    for k = 2:n
        % Prediction step
        x_pred(k) = x_filt(k-1);     % Predicted AoA
        P_pred(k) = P_filt(k-1) + Q; % Predicted covariance

        % Measurement update step
        K = P_pred(k) / (P_pred(k) + R);                       % Kalman gain
        x_filt(k) = x_pred(k) + K * (AoA_meas(k) - x_pred(k)); % Updated state
        P_filt(k) = (1 - K) * P_pred(k);                       % Updated covariance
    end

    % Backward pass (Rauch–Tung–Striebel (RTS) smoother)
    x_smooth = x_filt;   % Initialize smoothed state
    P_smooth = P_filt;   % Initialize smoothed covariance

    for k = n-1:-1:1
        C = P_filt(k) / (P_pred(k+1) + eps);                           % Smoother gain
        x_smooth(k) = x_filt(k) + C * (x_smooth(k+1) - x_pred(k+1));   % Smoothed state
        P_smooth(k) = P_filt(k) + C^2 * (P_smooth(k+1) - P_pred(k+1)); % Smoothed covariance
    end

    AoA_kalman = x_smooth; % Final smoothed AoA over time
end

%% Plot comparison of raw, moving average, and RTS-smoothed AoA
figure; 
hold on; 
grid on;
plot(t, AoA_meas, 'r.', 'DisplayName', 'Raw AoA');
plot(t, AoA_smooth, 'b--', 'DisplayName', 'Moving Average');
plot(t, AoA_kalman, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Kalman Filter');
xlabel('Time (s)');
ylabel('AoA (deg)');
legend;
title('AoA Tracking Comparison: Raw vs Moving Average vs Kalman');
ylim([-10 20]);

%% Physical parameters to be set

g = 9.81;                    % Acceleration due to gravity (m/s^2)
h_radar = 0.3;              % Height of the radar above the ground (m)
h_launch = 0.1;                % Launch height of the ball (m)
d_radar_to_launch = 1.3;     % Horizontal distance from radar to ball launch point (m)
d_launch_to_net = 4.6;       % Horizontal distance from launch point to net (m)
v0 = predictedVelocity;      % Initial ball speed estimated from Doppler (m/s)
theta_max = max(AoA_kalman); % Maximum observed elevation AoA at the net from Kalman filter (degrees)

%% Calculation of ball height at net
d_radar_to_net = d_radar_to_launch + d_launch_to_net;
h_net = h_radar + d_radar_to_net * tand(theta_max);

%% Solve for launch angle
launch = @(theta) - (g * d_launch_to_net^2) / (2 * v0^2 * cosd(theta)^2) + ...
                d_launch_to_net * tand(theta) + (h_launch - h_net);

options = optimset('Display','off'); % Suppress fsolve output for cleaner display
theta_launch = fsolve(launch, 45, options); % Initial guess 45 deg

%% Calculate approximated maximum height reached by sports ball
h_max = h_launch + (v0^2 * sind(theta_launch)^2) / (2*g);

%% Generate full trajectory of ball
% Compute flight time until height = 0 (ground)
% Solve for t when y(t) = 0 using quadratic form:
% y(t) = h_launch + v0*sin(theta)*t - 0.5*g*t^2 = 0

a = -0.5 * g;                   % Coefficient of t^2 term (gravity effect)
b = v0 * sind(theta_launch);    % Coefficient of t term (initial vertical velocity)
c = h_launch;                   % Constant term (launch height)

% Compute positive root of quadratic equation (since time cannot be negative)
t_flight = (-b - sqrt(b^2 - 4*a*c)) / (2*a);  

% Generate evenly spaced time points over the flight duration
t = linspace(0, t_flight, 200);

% Compute horizontal and vertical positions at each time point
x = v0 * cosd(theta_launch) * t;                             % Horizontal displacement (m)
y = h_launch + v0 * sind(theta_launch) * t - 0.5 * g * t.^2; % Vertical displacement (m)

%% Plot trajectory
figure; 
hold on; 
grid on;
plot(x, y, 'b-', 'LineWidth', 2);
xlabel('Horizontal distance [m]');
ylabel('Height [m]');
title('Full Ball Trajectory');

% Mark important points
plot(d_launch_to_net, h_net, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
plot(0, h_launch, 'go', 'MarkerSize', 8, 'MarkerFaceColor','g');
plot(x(y==max(y)), max(y), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');

% Add title with predicted initial velocity and launch angle
title(sprintf('Full Ball Trajectory (v0 = %.2f m/s, launch = %.2f°)', v0, theta_launch));

legend('Trajectory', 'Height at net', 'Launch point', 'Max height', 'Location','Best');

ylim([0, max(y)+0.5]);
xlim([0, max(x)+0.5]);

%% Print summary
fprintf('\n----- Ball Trajectory Summary -----\n');
fprintf('Launch speed (from Doppler): %.2f m/s\n', v0);
fprintf('Max Angle of Arrival: %.2f deg\n', theta_max);
fprintf('Launch angle: %.2f deg\n', theta_launch);
fprintf('Maximum height: %.2f m\n', h_max);
fprintf('Height at net: %.2f m\n', h_net);
fprintf('Total horizontal range until ground: %.2f m\n', max(x));
% fprintf('Total flight time: %.2f s\n', t_flight);
