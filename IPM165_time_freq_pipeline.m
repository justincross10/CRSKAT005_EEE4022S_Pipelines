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

%% Plot Raw Spectrogram
figure;
imagesc(time_axis, freq_axis_cropped, P_cropped);
axis xy;
colormap parula;
c = colorbar;
c.Label.String = 'Magnitude [dB]';
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Raw Spectrogram of Data');
dcm = datacursormode(gcf);
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');

%% Fixed clutter subtraction and median frame subtraction
clutter_map = P_lin_cropped(:,10);         % Fixed clutter in first frame
P_clean_lin = P_lin_cropped - clutter_map; % Subtract clutter power in first frame from power in other frames
P_clean_lin = max(P_clean_lin, 0);         % Clip negative values

bg_freq = median(P_clean_lin, 2); % Median across time for each frequency bin
P_medsub = P_clean_lin - bg_freq; % Background subtraction of median at each time frame from clutter subtracted bins
P_medsub = max(P_medsub, 0);      % Keep non-negative

% Convert to dB
P_clean_db = 10*log10(P_medsub + eps);

%% Spectrogram after clutter subtraction and median frame subtraction
figure;
imagesc(time_axis, freq_axis_cropped, P_clean_db);
axis xy;                     
colormap parula;            
c = colorbar;
c.Label.String = 'Magnitude [dB]';                   
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram After Clutter Subtraction and Median Background Subtraction');
dcm = datacursormode(gcf);   % Get data cursor object for current figure
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');  % Show coordinates in a small window

%% Percentile threshold for further clutter removal to enhance target

% Apply percentile threshold at each frequency bin to only maintain power values that
% lie above threshold at different time frames
perc_row = 95;   
P_thresh = zeros(size(P_medsub));
for r = 1:size(P_medsub,1)
    thr = prctile(P_medsub(r,:), perc_row);
    P_thresh(r, :) = P_medsub(r, :) .* (P_medsub(r,:) >= thr);
end
P_perc = 10*log10(P_thresh + eps);

%% Spectrogram after percentile theshold applied
figure;
imagesc(time_axis, freq_axis_cropped, P_perc);
axis xy; 
colormap parula; 
c = colorbar;
c.Label.String = 'Magnitude [dB]';
title('Spectrogram After Application of Percentile-Threshold');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

%% Apply 1D Median Filter to Suppress Small Spectral Artifacts
P_med = zeros(size(P_perc));
for c = 1:size(P_perc, 2)
    P_med(:,c) = medfilt1(P_perc(:,c), 3); % Use a 3-point median filter on each column
end

%% Plot spectrogram after median filtering
figure;
imagesc(time_axis, freq_axis_cropped, P_med);
axis xy; 
colormap parula; 
c = colorbar;
c.Label.String = 'Magnitude [dB]';
title('Spectrogram After Median Filter');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

%% Find high power detections remaining

P_med_lin = 10.^(P_med/10);            % Convert from dB to linear power
threshold = prctile(P_med_lin(:), 50); % Power value at 90th percentile
candidates = P_med_lin > threshold;    % Select power candidates above threshold

%% Highlight high power detections
figure; 
imagesc(time_axis, freq_axis_cropped, P_med); 
axis xy; 
colormap parula; 
c = colorbar;
c.Label.String = 'Magnitude [dB]';
title('High Power Detections');
hold on;

% Plot candidates as red dots
[high_freq, high_time] = find(candidates);  % Find candidates at time and frequency
plot(time_axis(high_time), freq_axis_cropped(high_freq), 'r.', 'MarkerSize', 10);

hold off;

%% Create an array of all candidate detections at time and frequency point

candidates_time = time_axis(high_time);           
candidates_freq = freq_axis_cropped(high_freq);   
candidates_list = [candidates_time(:), candidates_freq(:)];

%% Identify streaks of high power to further select potential ball candidates

time_tol = 0.002; % Maximum amount of time between candidates
freq_tol = 50;    % Maximum frequency jump between candidates

streaks = {};  % Cell array of streaks
used = false(size(candidates_list,1),1);

for i = 1:size(candidates_list,1)
    if used(i)
        continue;
    end
    
    % Start new streak
    streak_time = candidates_list(i,1);
    streak_freq = candidates_list(i,2);
    used(i) = true;
    
    % Try to grow streak by checking later points
    for j = i+1:size(candidates_list,1)
        if ~used(j)
            dt = abs(candidates_list(j,1) - streak_time(end)); % Check difference in time between candidates
            df = abs(candidates_list(j,2) - streak_freq(end)); % Check difference in frequency between candidates
            
            if dt <= time_tol && df <= freq_tol % Do differences obey tolerance value? If yes, add to streak array
                streak_time(end+1) = candidates_list(j,1);
                streak_freq(end+1) = candidates_list(j,2);
                used(j) = true;
            end
        end
    end
    
    streaks{end+1} = [streak_time(:), streak_freq(:)];
end

%% Plot detected streaks after streak analysis
figure;
imagesc(time_axis, freq_axis_cropped, P_med);
axis xy;
colormap parula;
c = colorbar;
c.Label.String = 'Magnitude [dB]';
hold on;
title('Spectrogram of All Streaks');
dcm = datacursormode(gcf);   % Get data cursor object for current figure
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');  % Show coordinates in a small window

for s = 1:numel(streaks)
    plot(streaks{s}(:,1), streaks{s}(:,2), '-o', 'Color', 'r', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
end
hold off;

%% Filter out very short streaks using a time-freq area threshold to detect ball streaks only

min_area = 5; % Time-freq area threshold
valid_streaks = {};

for s = 1:numel(streaks)
    t_vals = streaks{s}(:,1);
    f_vals = streaks{s}(:,2);

    t_start = min(t_vals);  % Start time of streak
    t_end = max(t_vals);    % End time of streak
    f_start = min(f_vals);  % Start frequency of streak
    f_end = max(f_vals);    % End frequency of streak
    area = (t_end - t_start) * (f_end - f_start);

    if area >= min_area
        valid_streaks{end+1} = streaks{s}; % Add streak if it meets threshold
    end
end



%% Plot valid streaks only
figure;
imagesc(time_axis, freq_axis_cropped, P_med);
axis xy;
colormap parula;
c = colorbar;
c.Label.String = 'Magnitude [dB]';
hold on;
title('Streaks after Time–Frequency Area Filtering');
dcm = datacursormode(gcf);  % Get data cursor object for current figure
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');  % Show coordinates in a small window

for s = 1:numel(valid_streaks)
    plot(valid_streaks{s}(:,1), valid_streaks{s}(:,2), '-o', 'Color', 'r', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
end
hold off;

%% Merge valid streaks that are within same energy region to isolate one object
merge_time_tol = 0.002;   % Time threshold to allow valid streak to merge
merge_freq_tol = 50;     % Freq threshold to allow valid streak to merge

% Sort streaks by start time
[~, idx] = sort(cellfun(@(s) min(s(:,1)), valid_streaks));
streaks_sorted = valid_streaks(idx);

merged_streaks = {};
i = 1;

while i <= numel(streaks_sorted)
    current = streaks_sorted{i};
    j = i + 1;

    while j <= numel(streaks_sorted)
        next_streak = streaks_sorted{j};
        time_gap = min(next_streak(:,1)) - max(current(:,1)); % Gap between current streak end and next start
        freq_overlap = min(max(current(:,2)), max(next_streak(:,2))) - max(min(current(:,2)), min(next_streak(:,2))); % Determine if there is overlap in frequency ranges of streaks

        % If frequency ranges overlap or are within threshold
        freq_met = (freq_overlap > 0) || (abs(freq_overlap) <= merge_freq_tol);

        if (time_gap <= merge_time_tol) && freq_met % Merge streaks if both conditions are met
            current = [current; next_streak];
            j = j + 1; % Check next streak
        else
            break;
        end
    end

    % Add merged streak
    merged_streaks{end+1} = current;

    % Move past all merged streaks
    i = j;
end

%% Plot merged streaks
figure;
imagesc(time_axis, freq_axis_cropped, P_med);
axis xy;
colormap parula;
c = colorbar;
c.Label.String = 'Magnitude [dB]';
hold on;
title('Streaks after merging');
dcm = datacursormode(gcf);   % Get data cursor object for current figure
set(dcm, 'Enable', 'on', 'DisplayStyle', 'window');  % Show coordinates in a small window
xlabel('Time (s)');
ylabel('Frequency (Hz)');
for s = 1:numel(merged_streaks)
    plot(merged_streaks{s}(:,1), merged_streaks{s}(:,2), 'o', 'Color','r', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
end
hold off;

%% Extract velocity information from valid streaks

fprintf('The possible ball detections are:\n');

for s = 1:numel(merged_streaks)
    streak_filtered = merged_streaks{s};

    if isempty(streak_filtered)
        continue;
    end

    % Extract times and Doppler frequencies
    t_streak = streak_filtered(:,1);
    f_doppler = streak_filtered(:,2);

    % Convert to velocity
    c = 3e8;              % Speed of light
    f_c = 24.125e9;       % Radar carrier frequency
    v_streak = (c * f_doppler) / (2 * f_c); % Convert Doppler frequency to velocity

    % Remove duplicate time entries while preserving data alignment
    [t_unique, uniq_idx] = unique(t_streak, 'stable');
    v_unique = v_streak(uniq_idx);

    % Extract maximum velocity for each valid streak
    [v_max, idx_max_raw] = max(v_unique);
    t_max = t_unique(idx_max_raw);

    fprintf('Streak %d: Time %.3f–%.3f s, Freq %.2f–%.2f Hz, v_max = %.2f m/s\n', ...
        s, min(t_streak), max(t_streak), min(f_doppler), max(f_doppler), v_max);
end
