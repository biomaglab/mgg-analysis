function a_pipelinepreprocMGG(inputFile, Fs)
% A_PIPELINEPREPROCMGG
% Multichannel pipeline for EGG/MGG .dat files with amplitude in picoTesla (pT).

    if nargin < 1 || isempty(inputFile)
        [file, folder] = uigetfile( ...
            {'*.dat;*.txt;*.csv;*.tsv', 'Numeric data'; '*.*', 'All files'}, ...
            'Select the .dat file');
        if isequal(file, 0)
            fprintf('No file selected.\n');
            return;
        end
        inputFile = fullfile(folder, file);
    end

    %% Configuration
    cfg.firstColumnIsTime = true;
    cfg.defaultFs = 20;
    cfg.nSigma = 6.0;                  % Conservative threshold: 6 MAD
    cfg.hampelWindowSeconds = 0.55;   % Target total Hampel window (~11 samples at 20 Hz)
    cfg.hampelHalfWindow = [];         % Empty = calculated automatically from hampelWindowSeconds
    cfg.defaultLowHz = 0.015;
    cfg.defaultHighHz = 0.2;
    cfg.butterOrder = 5;
    cfg.movementFrequencyCpm = 3.0;
    cfg.movementWindowCycles = 1.0;
    cfg.picoTeslaFactor = 400.0;       % Conversion factor to picoTesla (pT) [*400] 

    %% Loading
    M = read_numeric_compatible(inputFile);
    
    M = M(~all(isnan(M), 2), :);
    M = M(:, ~all(isnan(M), 1));

    if isempty(M) || size(M,1) < 2
        error('The file does not contain a valid numeric matrix.');
    end

    if cfg.firstColumnIsTime && size(M,2) >= 2
        t = M(:,1);
        X = M(:,2:end) * cfg.picoTeslaFactor; % Conversion to pT applied (*400)
        dt = diff(t);
        dt = dt(isfinite(dt) & dt > 0);

        if nargin < 2 || isempty(Fs)
            if isempty(dt)
                Fs = cfg.defaultFs;
                warning('Could not infer Fs from the time data. Using %.6g Hz.', Fs);
            else
                Fs = 1 / median(dt);
            end
        end
    else
        X = M * cfg.picoTeslaFactor; % Conversion to pT applied (*400)
        if nargin < 2 || isempty(Fs)
            Fs = cfg.defaultFs;
        end
        t = (0:size(X,1)-1)' / Fs;
    end

    if ~isfinite(Fs) || Fs <= 0
        error('Fs must be positive.');
    end

    fprintf('\nFile: %s\n', inputFile);
    fprintf('Samples: %d | channels: %d | Fs: %.8g Hz | duration: %.2f min\n', ...
        size(X,1), size(X,2), Fs, (t(end)-t(1))/60);

    %% Correct NaN/Inf only
    for ch = 1:size(X,2)
        good = isfinite(X(:,ch));
        if nnz(good) < 2
            error('Channel %d contains fewer than two valid samples.', ch);
        end
        if ~all(good)
            X(~good,ch) = interp1(find(good), X(good,ch), ...
                find(~good), 'linear', 'extrap');
        end
    end
    raw = X;

    %% Initial user choices
    cfg.outlierMethod = request_outlier_method();
    [cfg.lowHz, cfg.highHz] = request_user_band( ...
        Fs, cfg.defaultLowHz, cfg.defaultHighHz);

    %% 1. Robust median centering
    centered = bsxfun(@minus, raw, nanmedian_compatible(raw,1));
    movementHz = cfg.movementFrequencyCpm / 60;

    %% Iterative processing with inspection BEFORE saving
    accepted = false;
    attempt = 0;
    inspectionFigures = [];

    while ~accepted
        attempt = attempt + 1;

        if ~isempty(inspectionFigures)
            for kk = 1:numel(inspectionFigures)
                if ishandle(inspectionFigures(kk))
                    close(inspectionFigures(kk));
                end
            end
        end

        fprintf('\n============================================================\n');
        fprintf('PROCESSING / INSPECTION - ATTEMPT %d\n', attempt);
        fprintf('============================================================\n');
        fprintf('Outliers: %s | nSigma = %.3g', upper(cfg.outlierMethod), cfg.nSigma);
        if strcmpi(cfg.outlierMethod,'hampel')
            fprintf(' | target window = %.3g s', cfg.hampelWindowSeconds);
        end
        fprintf('\n');
        fprintf('Butterworth: %.6g-%.6g Hz (%.4g-%.4g cpm), order %d, zero-phase\n', ...
            cfg.lowHz, cfg.highHz, cfg.lowHz*60, cfg.highHz*60, cfg.butterOrder);

        %% 2. Hampel or MAD
        switch lower(cfg.outlierMethod)
            case 'hampel'
                if isempty(cfg.hampelHalfWindow)
                    nWin = max(3, round(cfg.hampelWindowSeconds * Fs));
                    if mod(nWin,2) == 0, nWin = nWin + 1; end
                    k = floor(nWin/2);
                else
                    k = max(1, round(cfg.hampelHalfWindow));
                    nWin = 2*k + 1;
                end
                fprintf('Hampel: window = %d samples (%.3f s), half-window k=%d, threshold=%.3g MAD\n', ...
                    nWin, nWin/Fs, k, cfg.nSigma);
                [despiked, hampelCount] = hampel_multichannel(centered, k, cfg.nSigma);
                fprintf('Hampel corrected %d of %d points (%.3f%%).\n', ...
                    hampelCount, numel(centered), 100*hampelCount/numel(centered));
            case 'mad'
                despiked = mad_filter_multichannel(centered, cfg.nSigma);
            otherwise
                error('Invalid outlier method.');
        end

        %% 3. User-defined Butterworth filter
        lowHz = cfg.lowHz;
        highHz = cfg.highHz;
        if ~(0 < lowHz && lowHz < highHz && highHz < Fs/2)
            error('Band %.6g-%.6g Hz is invalid for Fs=%.6g Hz.', lowHz, highHz, Fs);
        end

        [z,p,kFilt] = butter(cfg.butterOrder, [lowHz highHz]/(Fs/2), 'bandpass');
        sos = zp2sos(z,p,kFilt);
        bandpassed = zeros(size(despiked));
        for ch = 1:size(despiked,2)
            bandpassed(:,ch) = sosfiltfilt_compatible(sos, despiked(:,ch));
        end

        %% 4. LMMSE movement filter
        [movementClean, movementNoise] = movement_filter_lmmse( ...
            bandpassed, Fs, movementHz, cfg.movementWindowCycles);

        %% Visual inspection
        inspectionFigures = create_figures_inspection(t, raw, centered, despiked, ...
            bandpassed, movementClean, cfg, attempt);

        fprintf('\nThe figures are open for inspection.\n');
        fprintf('No output files have been saved yet.\n');

        choice = request_inspection_decision(cfg.outlierMethod);

        switch choice
            case 1
                accepted = true;
            case 2
                if strcmpi(cfg.outlierMethod,'hampel')
                    cfg.outlierMethod = 'mad';
                else
                    cfg.outlierMethod = 'hampel';
                end
                fprintf('\nMethod changed to %s. Reprocessing...\n', upper(cfg.outlierMethod));
            case 3
                [cfg.lowHz, cfg.highHz] = request_user_band(Fs, cfg.lowHz, cfg.highHz);
                fprintf('\nBand changed. Reprocessing...\n');
            case 4
                fprintf('\nProcessing canceled by the user.\n');
                if ~isempty(inspectionFigures)
                    for kk = 1:numel(inspectionFigures)
                        if ishandle(inspectionFigures(kk)), close(inspectionFigures(kk)); end
                    end
                end
                return;
        end
    end

    %% Output - FILE NAMING AND SAVING
    [folder, base, ~] = fileparts(inputFile);

    % Clean the name by removing unwanted substrings
    cleanBase = strrep(base, '_original_dec20Hz', '');
    cleanBase = strrep(cleanBase, '_pipeline_sem_ICA', '');

    % Check whether the name contains 'pos' using strfind (compatible with all MATLAB versions)
    if ~isempty(strfind(lower(cleanBase), 'pos'))
        suffix = '(pos)';
    else
        suffix = '(pre)';
    end

    outputDir = fullfile(folder, [cleanBase suffix]);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % SAVE FILES TO THE OUTPUT FOLDER
    save_stage(fullfile(outputDir,'01_raw.dat'), t, raw);                       % RAW
    save_stage(fullfile(outputDir,'02_median_centered.dat'), t, centered);       % MEAN CENTERED
    save_stage(fullfile(outputDir,'03_outlier_filtered.dat'), t, despiked);
    save_stage(fullfile(outputDir,'04_bandpass.dat'), t, bandpassed);            % BUTTERWORTH
    save_stage(fullfile(outputDir,'05_movement_filtered.dat'), t, movementClean);% MOVEMENT FILTERED
    save_stage(fullfile(outputDir,'05_estimated_movement_noise.dat'), t, movementNoise);

    save(fullfile(outputDir,'pipeline_complete.mat'), ...
        't','Fs','raw','centered','despiked','bandpassed', ...
        'movementClean','movementNoise','cfg');

    create_report(fullfile(outputDir,'report.txt'), inputFile, Fs, cfg, raw, movementClean);
    save_inspection_figures(inspectionFigures, outputDir);

    fprintf('\nPipeline completed successfully.\n');
    fprintf('Results saved to: %s\n', outputDir);
end

%% ================= HELPER FUNCTIONS =================

function method = request_outlier_method()
    fprintf('\nOutlier rejection:\n');
    fprintf('  1 - Local Hampel\n');
    fprintf('  2 - Global MAD\n');
    s = input('Select [Enter = 1]: ', 's');
    if isempty(strtrim(s)) || strcmp(strtrim(s),'1')
        method = 'hampel';
    elseif strcmp(strtrim(s),'2')
        method = 'mad';
    else
        error('Invalid option. Use 1 or 2.');
    end
end

function [lowHz, highHz] = request_user_band(Fs, defaultLowHz, defaultHighHz)
    fprintf('\nSet the Butterworth band-pass filter in Hz.\n');
    s1 = input(sprintf('Lower frequency [%.6g Hz]: ', defaultLowHz), 's');
    s2 = input(sprintf('Upper frequency [%.6g Hz]: ', defaultHighHz), 's');

    if isempty(strtrim(s1))
        lowHz = defaultLowHz;
    else
        lowHz = str2double(strrep(strtrim(s1), ',', '.'));
    end
    if isempty(strtrim(s2))
        highHz = defaultHighHz;
    else
        highHz = str2double(strrep(strtrim(s2), ',', '.'));
    end

    if ~isfinite(lowHz) || ~isfinite(highHz)
        error('Frequencies must be numeric.');
    end
    if ~(lowHz > 0 && lowHz < highHz && highHz < Fs/2)
        error('Invalid band.');
    end
end

function [Y, nReplaced] = hampel_multichannel(X, k, nSigma)
    [nSamples,~] = size(X);
    Y = X;
    nReplaced = 0;
    for i = 1:nSamples
        si = max(1,i-k);
        ei = min(nSamples,i+k);
        segment = X(si:ei,:);
        med = nanmedian_compatible(segment,1);
        robustSD = 1.4826 * nanmedian_compatible(abs(bsxfun(@minus,segment,med)),1);
        replace = abs(X(i,:)-med) > nSigma.*robustSD & robustSD > eps;
        Y(i,replace) = med(replace);
        nReplaced = nReplaced + nnz(replace);
    end
end

function Y = mad_filter_multichannel(X, nSigma)
    Y = X;
    med = nanmedian_compatible(X,1);
    dev = abs(bsxfun(@minus,X,med));
    robustSD = 1.4826 * nanmedian_compatible(dev,1);
    threshold = nSigma.*robustSD;
    for ch = 1:size(X,2)
        replace = dev(:,ch) > threshold(ch) & robustSD(ch) > eps;
        Y(replace,ch) = med(ch);
    end
end

function [cleaned, estimatedNoise] = movement_filter_lmmse(X, Fs, frequencyHz, windowCycles)
    if frequencyHz <= 0 || windowCycles <= 0
        error('Frequency and window length must be positive.');
    end

    [nSamples,~] = size(X);
    winSamples = max(3, round((windowCycles/frequencyHz)*Fs));
    halfWin = floor(winSamples/2);
    localMean = zeros(size(X));
    localVar = zeros(size(X));

    for i = 1:nSamples
        si = max(1,i-halfWin);
        ei = min(nSamples,i+halfWin);
        segment = X(si:ei,:);
        localMean(i,:) = mean(segment,1);
        localVar(i,:) = var(segment,1,1);
    end

    signalVariance = mean(localVar,1);
    numerator = max(bsxfun(@minus,localVar,signalVariance),0);
    denominator = bsxfun(@max,localVar,signalVariance);
    denominator = max(denominator,eps);
    gain = numerator ./ denominator;
    estimatedNoise = localMean + gain .* bsxfun(@minus,X,localMean);
    cleaned = bsxfun(@minus,X,estimatedNoise);
end

function save_stage(filename, t, X)
    dlmwrite(filename, [t X], 'delimiter','\t','precision','%.10g');
end

function create_report(filename, inputFile, Fs, cfg, raw, finalSignal)
    fid = fopen(filename,'w');
    if fid < 0, return; end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'EGG/MGG PIPELINE - NO ICA\n');
    fprintf(fid,'File: %s\n',inputFile);
    fprintf(fid,'Fs: %.10g Hz\n',Fs);
    fprintf(fid,'Channels: %d\n',size(raw,2));
    fprintf(fid,'Outlier method: %s\n',upper(cfg.outlierMethod));
    fprintf(fid,'Butterworth: %.6g-%.6g Hz (%.4g-%.4g cpm)\n', cfg.lowHz,cfg.highHz,cfg.lowHz*60,cfg.highHz*60);
    fprintf(fid,'Final signal: 05_movement_filtered.dat\n');
end

function choice = request_inspection_decision(method)
    fprintf('\nDECISION AFTER INSPECTING THE PLOTS:\n');
    fprintf('  1 - Accept and SAVE\n');
    fprintf('  2 - Switch outlier method and review again\n');
    fprintf('  3 - Change Butterworth band and review again\n');
    fprintf('  4 - Cancel without saving\n');

    while true
        s = input('Select [1-4]: ','s');
        choice = str2double(strtrim(s));
        if isfinite(choice) && any(choice == 1:4)
            return;
        end
        fprintf('Invalid option. Enter 1, 2, 3, or 4.\n');
    end
end

function figs = create_figures_inspection(t, raw, centered, despiked, bandpassed, movementClean, cfg, attempt)
    stages = { ...
        'Original signal', raw; ...
        'After median centering', centered; ...
        sprintf('After %s', upper(cfg.outlierMethod)), despiked; ...
        sprintf('Butterworth %.4g-%.4g Hz',cfg.lowHz,cfg.highHz), bandpassed; ...
        'After LMMSE', movementClean};

    figs = zeros(1,size(stages,1));
    for s = 1:size(stages,1)
        X = stages{s,2};
        nCh = size(X,2);
        figs(s) = figure_fullscreen_compat('Visible','on','Color','w', ...
            'Name',sprintf('Attempt %d - %s',attempt,stages{s,1}), ...
            'NumberTitle','off');

        for ch = 1:nCh
            ax = subplot(nCh,1,ch);
            plot(ax,t/60,X(:,ch),'LineWidth',0.75);
            ylabel(ax,sprintf('Ch %d (pT)',ch));
            grid(ax,'on');
            if ch == 1
                title(ax,sprintf('%s | %s | %.4g-%.4g Hz', ...
                    stages{s,1},upper(cfg.outlierMethod),cfg.lowHz,cfg.highHz));
            end
            if ch < nCh
                set(ax,'XTickLabel',[]);
            else
                xlabel(ax,'Time (min)');
            end
            xlim(ax,[t(1) t(end)]/60);
        end
        drawnow;
    end
end

function fig = figure_fullscreen_compat(varargin)
    fig = figure(varargin{:});
    try
        set(fig,'Units','normalized','OuterPosition',[0 0 1 1]);
    catch
        try
            set(fig,'Position',get(0,'ScreenSize'));
        catch
        end
    end
end

function save_inspection_figures(figs, outputDir)
    names = { ...
        'figure_01_raw.png', ...
        'figure_02_median_centered.png', ...
        'figure_03_outlier_filtered.png', ...
        'figure_04_bandpass.png', ...
        'figure_05_movement_filtered.png'};

    for k = 1:min(numel(figs), numel(names))
        if isgraphics(figs(k))
            print(figs(k), fullfile(outputDir, names{k}), '-dpng', '-r160');
        end
    end
end

function y = sosfiltfilt_compatible(sos, x)
    x = x(:);
    y = x;
    for sec = 1:size(sos,1)
        b = sos(sec,1:3);
        a = sos(sec,4:6);
        if abs(a(1)) < eps, error('Invalid SOS section.'); end
        b = b./a(1);
        a = a./a(1);
        y = filtfilt(b,a,y);
    end
end

function M = read_numeric_compatible(filename)
    if exist('readmatrix','file') == 2
        M = readmatrix(filename);
        return;
    end
    try
        M = dlmread(filename);
        return;
    catch
    end
    imported = importdata(filename);
    if isstruct(imported)
        if isfield(imported,'data')
            M = imported.data;
        else
            error('The file does not contain a numeric matrix.');
        end
    else
        M = imported;
    end
end

function y = nanmean_compatible(x, dim)
    if nargin < 2, dim = 1; end
    valid = ~isnan(x);
    count = sum(valid,dim);
    x(~valid) = 0;
    y = sum(x,dim)./max(count,1);
    y(count==0) = NaN;
end

function y = nanmedian_compatible(x, dim)
    if nargin < 2, dim = 1; end
    if dim ~= 1
        order = [dim,1:dim-1,dim+1:ndims(x)];
        x = permute(x,order);
    else
        order = 1:ndims(x);
    end
    sz = size(x);
    x2 = reshape(x,sz(1),[]);
    y2 = NaN(1,size(x2,2));
    for kk = 1:size(x2,2)
        v = x2(:,kk);
        v = v(~isnan(v));
        if ~isempty(v), y2(kk) = median(v); end
    end
    outsz = sz; outsz(1)=1;
    y = reshape(y2,outsz);
    if dim ~= 1
        y = ipermute(y,order);
    end
end