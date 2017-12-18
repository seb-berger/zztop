function analyse()
%ANALYSE reproduces the results presented in the article [1].
%   To carry out the analyses, please obtain the EEG datasets 'CAP Sleep
%   Database' and 'CHB-MIT Scalp EEG Database', which are freely available
%   from PhysioNet (www.physionet.org).
%
%   By default, the function expects the corresponding EDF files in the
%   subdirectories 'cap' and 'chb_mit'. This can be changed by modifying
%   the variables 'cap_path' and 'mit_path' in the code below.
%
%   Additional settings can also be tweaked: please see the 'Constants'
%   section in the code of this function.
%
%   The function generates a considerable amount of intermediate results
%   and stores them in the 'cap' and 'chb_mit' subdirectories. Please
%   arrange for approximately 50 GByte of free storage space.
%
%   Results are produced in the order in which they appear in the
%   manuscript. They are displayed on screen, and additionally, LaTeX code
%   is generated. Thus, the exact origins of all data reported in [1]
%   should be reproducible. Notice, however, that different environments
%   (GNU Octave, MATLAB, etc.) and software versions may produce slightly
%   different results. Data presented in [1] were generated using
%   GNU Octave, version 4.2.1.
%
%   Thank you for your interest in this publication. In particular,
%   your looking into the details of the code is highly appreciated!
%
%   Sebastian Berger, 2017.
%
%   [1] Berger, S. Permutation Entropy: Too Complex a Measure For EEG Time
%       Series? Entropy 2017, xx, xxxx-xxxx.

%   Copyright (c) 2017, Sebastian Berger.
%
%   Klinikum rechts der Isar der
%   Technischen Universitaet Muenchen
%   Munich, Germany
%
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are
%   met:
%       * Redistributions of source code must retain the above copyright
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright
%         notice, this list of conditions and the following disclaimer in
%         the documentation and/or other materials provided with the
%         distribution.
%       * Neither the names of the copyright holders nor the names of its
%         contributors may be used to endorse or promote products derived
%         from this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
%   THE KLINIKUM RECHTS DER ISAR BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
%   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
%   THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.


%%% Platform dependent preamble %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isoctave()
    % Load packages and enable unbuffered screen output for GNU Octave.
    pkg('load', 'signal');
    pkg('load', 'statistics');
    page_screen_output(false, 'local');
else
    % Create function aliases for MATLAB.
    princomp = @(x) pca(x);
    spearman = @(x, y) corr(x, y, 'type', 'Spearman');
end


%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Whether to keep existing EEG processing data
overwrite_old = false;

% Path to CAP Sleep Database
cap_path = 'cap';

% Path to CHB_MIT Scalp EEG Database
mit_path = 'chb_mit';

% Path to store results to
dest_path = fullfile('.', 'data');

% Channels to be used for the CAP Sleep Database
cap_labels = {'C3', 'O2',  'FP1',  'O2A1', 'F3-C3',  'P4-O2', ...
              'C4', 'P3',  'Fp2', 'C3-A2', 'F4-C4',  'T3-T5', ...
              'F3', 'P4', 'C3A2', 'C3-P3', 'F7-T3',  'T4-T6', ...
              'F4', 'T3', 'C4A1', 'C4-A1', 'F8-T4', 'FP1-F3', ...
              'F7', 'T4', 'F3A2', 'C4-P4', 'O1-A2', 'FP2-F4,' ...
              'F8', 'T5', 'F4A1', 'F1-F3', 'O2-A1', 'Fp2-F4', ...
              'O1', 'T6', 'O1A2', 'F2-F4', 'P3-O1'};

% Channels to be used for the CHB-MIT Scalp EEG Database
mit_labels = {'C2',       'C2-CS2',  'C3',      'C3-CS2',  ...
              'C3-P3',    'C4',      'C4-CS2',  'C4-P4',   ...
              'C6',       'C6-CS2',  'CP1-Ref', 'CP2',     ...
              'CP2-CS2',  'CP2-Ref', 'CP4',     'CP4-CS2', ...
              'CP5-Ref',  'CP6',     'CP6-CS2', 'CP6-Ref', ...
              'CZ',       'CZ-CS2',  'CZ-PZ',   'F3',      ...
              'F3-C3',    'F3-CS2',  'F4',      'F4-C4',   ...
              'F4-CS2',   'F7',      'F7-CS2',  'F7-T7',   ...
              'F8',       'F8-CS2',  'F8-T8',   'FC1-Ref', ...
              'FC2-Ref',  'FC5-Ref', 'FC6-Ref', 'FP1',     ...
              'FP1-CS2',  'FP1-F3',  'FP1-F7',  'FP2',     ...
              'FP2-CS2',  'FP2-F4',  'FP2-F8',  'FT10-T8', ...
              'FT9-FT10', 'FZ',      'FZ-CS2',  'FZ-CZ',   ...
              'O1-CS2',   'O2',      'O2-CS2',  'P3',      ...
              'P3-CS2',   'P3-O1',   'P4',      'P4-CS2',  ...
              'P4-O2',    'P7',      'P7-CS2',  'P7-O1',   ...
              'P7-T7',    'P8',      'P8-CS2',  'P8-O2',   ...
              'PZ',       'PZ-CS2',  'PZ-OZ',   'T7',      ...
              'T7-CS2',   'T7-FT9',  'T7-P7',   'T8',      ...
              'T8-CS2',   'T8-P8'};

% Kernel density bandwith for PeEn
pe_density_sigma  = 5e-3;

% Kernel density bandwith for PCA scores
pca_density_sigma = 1e-2;


%%% EEG Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print section header
fprintf('%s\n', '=' * ones(1, 80));
fprintf('%s', ' ' * ones(1, 24));
fprintf('Step 1 of 2: Process EEG Data.\n');
fprintf('%s\n\n', '=' * ones(1, 80));

% Process the CAP Sleep Database
fprintf('Processing the CAP Sleep Database...\n');
if isdir(cap_path)
    files = findfiles('*.edf', cap_path);

    if ~isempty(files)
        for f = 1:numel(files)
            [p, n] = fileparts(files{f});
            res_filename = fullfile(p, [n, '.res.mat']);

            % Check if results file already exists
            if ~overwrite_old && exist(res_filename, 'file')
                continue;
            end

            % Analyse EDF file
            results = process_edf_file(files{f}, cap_labels, 200, 4000);

            % Store results
            save(res_filename, 'results', '-v6');
        end
    else
        fprintf('No EDF files found in directory ''%s''.\n', cap_path);
    end
else
    fprintf('Directory ''%s'' does not exist.\n', cap_path);
end

% Process the CHB-MIT Scalp EEG Database
fprintf('Processing the CHB-MIT Scalp EEG Database...\n');

if isdir(mit_path)
    files = findfiles('*.edf', mit_path);

    if ~isempty(files)
        for f = 1:numel(files)
            [p, n] = fileparts(files{f});
            res_filename = fullfile(p, [n, '.res.mat']);

            % Check if results file already exists
            if ~overwrite_old && exist(res_filename, 'file')
                continue;
            end

            % Analyse EDF file
            results = process_edf_file(files{f}, mit_labels, 200, 4000);
            if isempty(results)
                a = 1;
            end

            % Store results
            save(res_filename, 'results', '-v6');
        end
    else
        fprintf('No EDF files found in directory ''%s''.\n', mit_path);
    end
else
    fprintf('Directory ''%s'' does not exist.\n', mit_path);
end


%%% Data analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print section header
fprintf('\n%s\n', '=' * ones(1, 80));
fprintf('%s', ' ' * ones(1, 24));
fprintf('Step 2 of 2: Analyse results.\n');
fprintf('%s\n\n', '=' * ones(1, 80));

% Prepare subdirectory
if ~isdir(fullfile('.', dest_path))
    mkdir(fullfile('.', dest_path));
end

% Load results data
fprintf('Loading data...\n');

cap_files = findfiles('*.res.mat', cap_path);
mit_files = findfiles('*.res.mat', mit_path);

cap_res = load_results(cap_files, 3);
mit_res = load_results(mit_files, 3);

% Concatenate results
cap_peaks    = [cap_res.peaks];    % Number of peaks
cap_pe       = [cap_res.pe_ent];   % PeEn
cap_peak_ent = [cap_res.peak_ent]; % Entropy of peaks
balances     = [cap_res.bal];      % Balance coefficients
mit_peaks    = [mit_res.peaks];    % Number of peaks
mit_pe       = [mit_res.pe_ent];   % PeEn

% Calculate peak probabilities
cap_peak_p   = cap_peaks / 3998;
mit_peak_p   = mit_peaks / 3998;

% Get number of EEG epochs
cap_num      = size(cap_pe, 2);
mit_num      = size(mit_pe, 2);

% Calculate PeEn kernel density for CAP data
fprintf('Calculating PeEn density...\n');
pe_den_idx = linspace(0, 1 + 3 * pe_density_sigma, 200);
[pe_den, pe_den_idx] = kden(cap_pe, pe_density_sigma, pe_den_idx);

% Create a new figure
screen_size = get(0, 'ScreenSize');
screen_size = screen_size(3:4);

% Make sure figure fits screen
dims = [min(1280, screen_size(1)), min(1024, screen_size(2))];
pos = (screen_size - dims) / 2;
figure('Position', [pos, dims]);

% Visualise PeEn density
subplot(3, 4, 1:2);
plot(pe_den_idx, pe_den);
title('PeEn');
xlabel('PeEn / log6');
ylabel('Density');
drawnow();

% Store PeEn kernel density in data file
fid = fopen(fullfile(dest_path, 'peen_density.dat'), 'w');
res = [pe_den_idx; pe_den];
fprintf(fid, '%g %g\n', res(:));
fclose(fid);

% Calculate PCA on the balance coefficients
fprintf('Calculating PCA...\n');
[weights, scores, eigenvals] = princomp(balances.');

% Sign is implementation-dependent: flip if necessary!
if weights(1, 1) > 0
    weights = -1 * weights;
    scores  = -1 * scores;
end

% Calculate principle component kernel densities
for n = 1:4
    fprintf('Calculating PCA density %d...\n', n);
    if n == 1
        limits = linspace(-0.5, 1, 200);
    else
        limits = linspace(-0.2, 0.2, 200);
    end

    [den, idx] = kden(scores(:, n), pca_density_sigma, limits);

    % Visualise kernel density
    subplot_idx = [5, 6, 9, 10];
    subplot(3, 4, subplot_idx(n));
    plot(idx, den);
    title(['z', num2str(n)]);
    xlabel('Score');
    ylabel('Density');
    drawnow();

    % Store result in data file
    filename = fullfile(dest_path, ['pca_density_', num2str(n), '.dat']);
    fid = fopen(filename, 'w');
    res = [idx; den];
    fprintf(fid, '%g %g\n', res(:));
    fclose(fid);
end

fprintf('Calculating PCA Correlations...\n');

% Calculate explained variances
explained = eigenvals / sum(eigenvals);

% Calculate Spearman correlations
fprintf('============================================\n');
fprintf(' Component  Eigenvalue  Explained  Spearman \n');
fprintf('============================================\n');

% Create a TeX file
fid = fopen(fullfile(dest_path, 'spearman_balance_corr_data.tex'), 'w');

for n = 1:4
    % Correlation between PeEn and principle component n
    spear = spearman(scores(:, n), cap_pe.');

    % Print results to on-screen table
    fprintf('%10s  ', ['z', num2str(n)]);
    fprintf('%10.1d  ', eigenvals(n));
    fprintf('%8.1f%%  ', 100 * explained(n));
    fprintf('%8.5f\n', spear);

    % Store results in TeX file
    fprintf(fid, '$\\mathbf{z}_%d$ & ', n);
    fprintf(fid, '$\\scinum{%.3e}{2}$ & ', eigenvals(n));
    fprintf(fid, '$%.1f \\%%$ & ', 100 * explained(n));

    if n == 1
        fprintf(fid, '$%.5g$ \\\\\n', spear);
    else
        fprintf(fid, '$%.1g$ \\\\\n', spear);
    end
end
fclose(fid);
fprintf('============================================\n\n');

% Calculate statistics of balance coefficients
fprintf('Calculating mean, median and mode of balance coefficients...\n');
patterns = sortrows(perms(1:3));
pairs = nchoosek(1:6, 2);

means = mean(balances, 2);
medians = median(balances, 2);
modes = mode(balances, 2);

% Create TeX file
fid = fopen(fullfile(dest_path, 'first_prin_comp_data.tex'), 'w');

fprintf('========================================================\n');
fprintf(' Index  Balance       Weight       Mean     Med    Mode \n');
fprintf('========================================================\n');

% First print peak-to-non-peak balances, then the remaining patterns
for n = [1, 2, 3, 4, 9, 12, 14, 15, 5, 6, 7, 8, 10, 11, 13]
    % Insert a separator before index 5
    if n == 5
        fprintf('----------------------------');
        fprintf('----------------------------\n');
        fprintf(fid, '\\midrule\n');
    end

    % Get string representations of the patterns
    pat1 = sprintf('%d', patterns(pairs(n, 1), :));
    pat2 = sprintf('%d', patterns(pairs(n, 2), :));

    % Print table row
    fprintf('%6d  '  , n);
    fprintf('%7s  '  , [pat1, '/', pat2]);
    fprintf('%6.3f  ', weights(n, 1));
    fprintf('%6.3f  ', weights(n, 1) / weights(1, 1));
    fprintf('%6.4f  ', means(n));
    fprintf('%6.4f  ', medians(n));
    fprintf('%6.4f\n', modes(n));

    % Store result in TeX file
    fprintf(fid, '%%\n');
    fprintf(fid, '%d &\n', n);
    fprintf(fid, '$\\balcoeff{%s}{%s}$ &\n', pat1, pat2);
    fprintf(fid, '$%.3f$ &\n', weights(n, 1));
    fprintf(fid, '$%.3f$ &\n', weights(n, 1) / weights(1, 1));
    fprintf(fid, '$%.4f$ &\n', means(n));
    fprintf(fid, '$%.4f$ &\n', medians(n));
    fprintf(fid, '$%.4f$ \\\\\n', modes(n));
end
fclose(fid);
fprintf('========================================================\n\n');

fprintf('Calculating Pearson correlation...\n');

% Correlate PeEn and entropy of peaks
pearson = corr(cap_pe.', cap_peak_ent.');
rel_err = nanmean(abs(cap_pe - cap_peak_ent) ./ cap_pe);

% Print results
fprintf('Pearson correlation between PeEn and EoP: %.6f\n', pearson);
fprintf('Relative error: %.2g%%.\n', 100 * rel_err);

fprintf('Calculating Spearman correlations...\n');

% Split CAP data into subsets
cap_pe_lo     = cap_pe(cap_peak_p <= 2/3);
cap_pe_hi     = cap_pe(cap_peak_p >= 2/3);
cap_peak_p_lo = cap_peak_p(cap_peak_p <= 2/3);
cap_peak_p_hi = cap_peak_p(cap_peak_p >= 2/3);

% Split CHB_MIT data into subsets
mit_pe_lo     = mit_pe(mit_peak_p <= 2/3);
mit_pe_hi     = mit_pe(mit_peak_p >= 2/3);
mit_peak_p_lo = mit_peak_p(mit_peak_p <= 2/3);
mit_peak_p_hi = mit_peak_p(mit_peak_p >= 2/3);

% Calculate Spearman correlations between PeEn and peak probabilities
cap_spear    = spearman(cap_pe.', cap_peak_p.');
cap_spear_lo = spearman(cap_pe_lo.', cap_peak_p_lo.');
cap_spear_hi = spearman(cap_pe_hi.', cap_peak_p_hi.');
mit_spear    = spearman(mit_pe.', mit_peak_p.');
mit_spear_lo = spearman(mit_pe_lo.', mit_peak_p_lo.');
mit_spear_hi = spearman(mit_pe_hi.', mit_peak_p_hi.');

% Create on-screen table
fprintf('================================================================\n');
fprintf('             Number of EEG Epochs        Spearman Correlation   \n');
fprintf('Database  -------------------------  ---------------------------\n');
fprintf('              all       lo       hi      all       lo        hi \n');
fprintf('================================================================\n');
fprintf(' CAP      ');
fprintf('%7.2g  ', numel(cap_pe));
fprintf('%7.2g  ', numel(cap_pe_lo));
fprintf('%7.2g  ', numel(cap_pe_hi));
fprintf('%7.5f  ', cap_spear);
fprintf('%7.5f  ', cap_spear_lo);
fprintf('%7.5f\n', cap_spear_hi);

fprintf(' CAP-MIT  ');
fprintf('%7.2g  ', numel(mit_pe));
fprintf('%7.2g  ', numel(mit_pe_lo));
fprintf('%7.2g  ', numel(mit_pe_hi));
fprintf('%7.5f  ', mit_spear);
fprintf('%7.5f  ', mit_spear_lo);
fprintf('%7.5f\n', mit_spear_hi);

fprintf('================================================================\n\n');

% Store results in TeX file
fid = fopen(fullfile(dest_path, 'correlations_data.tex'), 'w');
fprintf(fid, 'CAP &\n');
fprintf(fid, '$\\scinum{%g}{2}$ &\n', numel(cap_pe));
fprintf(fid, '$\\scinum{%g}{2}$ &\n', numel(cap_pe_lo));
fprintf(fid, '$\\scinum{%g}{2}$ &\n', numel(cap_pe_hi));
fprintf(fid, '$%.5g$ &\n', cap_spear);
fprintf(fid, '$%.5g$ &\n', cap_spear_lo);
fprintf(fid, '$%.3g$ \\\\\n', cap_spear_hi);

fprintf(fid, 'CHB-MIT &\n');
fprintf(fid, '$\\scinum{%g}{2}$ &\n', numel(mit_pe));
fprintf(fid, '$\\scinum{%g}{2}$ &\n', numel(mit_pe_lo));
fprintf(fid, '$\\scinum{%g}{2}$ &\n', numel(mit_pe_hi));
fprintf(fid, '$%.5g$ &\n', mit_spear);
fprintf(fid, '$%.5g$ &\n', mit_spear_lo);
fprintf(fid, '$%.3g$ \\\\\n', mit_spear_hi);
fclose(fid);

% Calculate accumulated explained variations
explained1 = round(explained(1) * 100 - 0.5);
explained4 = round(sum(explained(1:4)) * 1000 - 0.5) / 10;

% Store remaining constants in a common TeX file
fid = fopen(fullfile(dest_path, 'constants.tex'), 'w');
fprintf(fid, '\\newcommand{\\ncap}{%d}\n', cap_num);
fprintf(fid, '\\newcommand{\\nmit}{%d}\n', mit_num);
fprintf(fid, '\\newcommand{\\peendensitysigma}{%g}\n', pe_density_sigma);
fprintf(fid, '\\newcommand{\\pcaexplained}{%g}\n', explained4 / 100);
fprintf(fid, '\\newcommand{\\pcaexplainedfirstpercent}{%g}\n', explained1);
fprintf(fid, '\\newcommand{\\pcaexplainedpercent}{%g}\n', explained4);
fprintf(fid, '\\newcommand{\\pearsoncorr}{%.5g}\n', pearson);
fprintf(fid, '\\newcommand{\\relerr}{%.1g}\n', rel_err * 100);
fclose(fid);

% Prepare the preliminary results plots for higher pattern orders
for ord = 3:5
    fprintf('Calculating percentiles for order %d...\n', ord);

    % Data for order 3 is still in memory!
    if ord > 3
        cap_res = load_results(cap_files, ord);
        mit_res = load_results(mit_files, ord);
    end

    % Concatenate results from both databases
    peaks = [cap_res.peaks,  mit_res.peaks];
    pe    = [cap_res.pe_ent, mit_res.pe_ent];

    % Obtain 101 percentiles, omit bins that are less than 100
    [res, idx] = percentiles(peaks, pe, 101, 100);

    % Scale peak count to peak probability
    idx = idx / 3998;

    % Select 5th, 25th, 50th, 75th and 95th percentile
    res = res(:, [6, 26, 51, 76, 96]);

    % Apply smoothing
    res = conv2(res, ones(10, 1) / 10, 'same');

    % Omit edges (where smoothing kernel did not fully overlap)
    idx = idx(5:end-5);
    res = res(5:end-5, :);

    % Plot percentile bands
    subplot_idx = [3, 4; 7, 8; 11, 12];
    subplot(3, 4, subplot_idx(ord - 2, :));
    palette = get(gca(), 'ColorOrder');
    color = palette(1, :);

    plot(idx, res, 'Color', color);
    title(['Order m = ', num2str(ord)]);
    xlabel('Peak probability');
    ylabel('PeEn');
    drawnow();

    % Store data
    dlmwrite(fullfile(dest_path, sprintf('percentiles_m%d.dat', ord)), ...
        [idx, res], ' ');
end


%%% Auxiliary function: process_edf_file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = process_edf_file(filename, channels, fs, wnd)

res = [];

% Print header
[~, n, e] = fileparts(filename);
fprintf('\n%s\n', '=' * ones(1, 80));
fprintf(' %s\n', [n, e]);
fprintf('%s\n', '=' * ones(1, 80));

% Try to load the file
fprintf(' Loading EDF file...                    ');
try
    edf = import_edf(filename, false);
catch err
    fprintf('Error: %s\n\n', err.message);
    return
end
fprintf('Done.\n');

% Delete unneeded channels
mask = ismember({edf.chan.label}, channels);
edf.chan = edf.chan(mask);
edf.num_ch = sum(mask);

% Gather the set of sample rates
rates = [edf.chan.num_smp] / edf.dur_rcd;

% Resample all channels to common sample rate fs
fprintf(' Resampling to %d Hz...                ', fs);
fcn = @(r) par_resample([edf.chan(rates == r).samples], fs, r);
eeg = arrayfun(fcn, unique(rates), 'UniformOutput', false);
eeg = [eeg{:}];
fprintf('Done.\n');

% Count peaks
fprintf(' Counting peaks...                      ');
peak_count = abs(diff(diff(eeg) >= 0));
peak_count = peak_count(1:floor(end / wnd)*wnd, :);
peak_count = reshape(peak_count, wnd, []);
peak_count = sum(peak_count(1:end-2, :));
fprintf('Done.\n');

for ord = 3:5
    % Encode EEG into ordinal patterns
    fprintf(' Extracting patterns of order m = %d...  ', ord);
    patterns = symbolise(eeg, ord, 1, true);

    % Split patterns into blocks that correspond to 'wnd' EEG samples
    patterns = patterns(1:(wnd * floor(end / wnd)), :);
    patterns = reshape(patterns, wnd, []);
    patterns = patterns(1:(end - ord + 1), :);

    % Remove any blocks that contain ties
    tie_mask = any(patterns == -1);
    patterns = patterns(:,   ~tie_mask);
    peaks    = peak_count(:, ~tie_mask);

    % Store peak count
    tmp = [];
    tmp.peaks = peaks;
    fprintf('Done.\n');

    % Count the number of occurrences for each pattern type
    fprintf(' Estimating pattern probabilities...    ');
    probs = histc(patterns, 1:factorial(ord)) / (wnd - ord + 1);
    fprintf('Done.\n');

    % Compute PeEn
    fprintf(' Calculating permutation entropy...     ');
    tmp.pe_ent = -nansum(probs .* log2(probs)) / log2(factorial(ord));
    fprintf('Done.\n');

    if ord == 3
        % Compute entropy of peaks
        fprintf(' Calculating entropy of peaks...        ');
        peak = sum(probs(2:5, :));
        ent = peak .* log2(peak) + (1 - peak) .* log2(1 - peak);
        ent(isnan(ent)) = 0;
        tmp.peak_ent = (peak + 1 - ent) / log2(6);
        fprintf('Done.\n');

        % Compute balance coefficients
        fprintf(' Calculating balance coefficients...    ');
        pairs = nchoosek(1:6, 2);

        tmp.bal = probs(pairs(:, 1), :) ./ ...
            (probs(pairs(:, 1), :) + probs(pairs(:, 2), :));

        tmp.bal(isnan(tmp.bal)) = 0.5;
        fprintf('Done.\n');
    end

    % Save number of ties
    tmp.num_ties = sum(tie_mask);

    % Store in struct to be returned
    res.(sprintf('m%d', ord)) = tmp;
end


%%% Auxiliary function: load_results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = load_results(files, order)

for f = 1:numel(files)
    tmp = load(files{f});
    tmp = tmp.results;
    if ~isempty(tmp)
        break;
    end
end

num_omitted = f - 1;

tmp = tmp.(sprintf('m%d', order));
fields = fieldnames(tmp).';
fields = [fields; cell(size(fields))];
res = struct(fields{:});

res(f) = tmp;
for f = f+1:numel(files)
    tmp = load(files{f});
    tmp = tmp.results;
    if isempty(tmp)
        num_omitted = num_omitted + 1;
        continue;
    end

    res(f) = tmp.(sprintf('m%d', order));
end
