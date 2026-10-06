function c_quantify_metrics(preFile, postFile, userMinFilterCpm, userMaxFilterCpm)
% C_QUANTIFY_METRICS
% Processes and compares recordings of up to 10 minutes (pre- and postprandial),
% with user-defined filter bands and GA-RI calculation,
% with optional BMI-based adjustment of amplitude and GA-RI.
%
% Linear detrending is applied to the global signal to remove baseline offsets,
% drift and low-frequency artifacts from the PSD.

    if nargin < 2 || isempty(preFile) || isempty(postFile)
        [selectedPreFile, preFolder] = uigetfile({'*.dat;*.txt;*.csv','Data files'}, 'Select the PRE-PRANDIAL file');
        if isequal(selectedPreFile,0), return; end
        preFile = fullfile(preFolder, selectedPreFile);
        
        [selectedPostFile, postFolder] = uigetfile({'*.dat;*.txt;*.csv','Data files'}, 'Select the POST-PRANDIAL file');
        if isequal(selectedPostFile,0), return; end
        postFile = fullfile(postFolder, selectedPostFile);
    end

    %% 1. Global settings and user input (MGG filter + BMI correction)
    cfg.samplingRate = 20.0;             % Default sampling frequency (Hz)
    cfg.nfft = 32768;
    cfg.picoTeslaFactor = 1.0;   % Unit factor: the input signal is already in pT
    cfg.figureDpi = 300;       % Figure resolution

    % Request the band-pass filter frequencies and BMI value
    if nargin < 4
        prompt = {'What was the minimum frequency of the band-pass filter used in the MGG pipeline? (cpm):', ...
                  'What was the maximum frequency of the band-pass filter used in the MGG pipeline? (cpm):', ...
                  'Enter the patient BMI (enter 0 or leave blank to NOT apply adjustment):'};
        dialogTitle = 'Input Parameters and Anatomical Adjustment';
        dims = [1 75];
        defaultInput = {'1.0', '6.0', '0'}; % Recommended default remains 1.0 - 6.0 cpm
        userResponse = inputdlg(prompt, dialogTitle, dims, defaultInput);
        if isempty(userResponse), return; end
        cfg.filterMinCpm = str2double(userResponse{1});
        cfg.filterMaxCpm = str2double(userResponse{2});
        bmiValue       = str2double(userResponse{3});
    else
        cfg.filterMinCpm = userMinFilterCpm;
        cfg.filterMaxCpm = userMaxFilterCpm;
        
        prompt = {'Enter the patient BMI (enter 0 to NOT apply adjustment):'};
        defaultInput = {'0'};
        userResponse = inputdlg(prompt, 'BMI adjustment', [1 60], defaultInput);
        if isempty(userResponse)
            bmiValue = 0;
        else
            bmiValue = str2double(userResponse{1});
        end
    end

    % Process BMI input
    if isnan(bmiValue) || bmiValue <= 0
        cfg.bmi = 0;
        cfg.hasBmi = false;
    else
        cfg.bmi = bmiValue;
        cfg.hasBmi = true;
    end

    % Convert the frequency limits from cpm to Hz
    cfg.fMinHz = cfg.filterMinCpm / 60.0;
    cfg.fMaxHz = cfg.filterMaxCpm / 60.0;

    %% 2. Process the signals (pre- and postprandial)
    fprintf('\n---> Processing PRE-PRANDIAL data (band: %.1f to %.1f cpm)...\n', cfg.filterMinCpm, cfg.filterMaxCpm);
    preResults = process_record(preFile, cfg);
    
    fprintf('\n---> Processing POST-PRANDIAL data (band: %.1f to %.1f cpm)...\n', cfg.filterMinCpm, cfg.filterMaxCpm);
    postResults = process_record(postFile, cfg);

    %% 3. Apply BMI correction (if provided)
    if cfg.hasBmi
        bmiFactor = exp(-0.05 * (cfg.bmi - 23.8));
        
        preResults.meanAmplitude = preResults.meanAmplitude / bmiFactor;
        preResults.amplitudeP90   = preResults.amplitudeP90 / bmiFactor;
        postResults.meanAmplitude = postResults.meanAmplitude / bmiFactor;
        postResults.amplitudeP90   = postResults.amplitudeP90 / bmiFactor;
        
        preResults.gari = min(1.0, max(0.0, preResults.gari / bmiFactor));
        postResults.gari = min(1.0, max(0.0, postResults.gari / bmiFactor));
        
        fprintf('[INFO] BMI correction (%.1f kg/m^2) applied successfully.\n', cfg.bmi);
    else
        fprintf('[INFO] Processing completed with raw values (no BMI adjustment).\n');
    end

    %% 4. Calculate combined metrics (BSGM / EGG)
    amplitudeRatio = postResults.meanAmplitude / preResults.meanAmplitude;
    fedFastedAmplitudeRatio    = postResults.amplitudeP90 / preResults.meanAmplitude;

    %% 5. Display results in the command window
    fprintf('\n======================================================\n');
    fprintf('           PRE vs POST ANALYSIS RESULTS               \n');
    fprintf('======================================================\n');
    fprintf('MGG filter band     : %.2f - %.2f cpm\n', cfg.filterMinCpm, cfg.filterMaxCpm);
    if cfg.hasBmi
        fprintf('BMI adjustment      : YES (BMI = %.1f kg/m^2)\n', cfg.bmi);
    else
        fprintf('BMI adjustment      : NO (raw values)\n');
    end
    fprintf('Pre-prandial GA-RI  : %.4f (normal >= 0.25)\n', preResults.gari);
    fprintf('Post-prandial GA-RI : %.4f (normal >= 0.25)\n', postResults.gari);
    fprintf('Pre mean amplitude  : %.4f pT\n', preResults.meanAmplitude);
    fprintf('Post mean amplitude : %.4f pT\n', postResults.meanAmplitude);
    fprintf('Post P90 amplitude  : %.4f pT\n', postResults.amplitudeP90);
    fprintf('Amplitude Ratio     : %.4f\n', amplitudeRatio);
    fprintf('ff-AR (Fed:Fasted)  : %.4f\n', fedFastedAmplitudeRatio);
    fprintf('======================================================\n');

    %% 6. Plot the results
    fig = plot_comparison(preResults, postResults, amplitudeRatio, fedFastedAmplitudeRatio, cfg);

    %% 7. Save the figure and data table
    [folderPre, namePre, ~] = fileparts(preFile);
    outDir = fullfile(folderPre, [namePre '_comparativo_bsgm']);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    figFilePng = fullfile(outDir, 'comparativo_pre_pos_bsgm.png');
    print(fig, figFilePng, '-dpng', sprintf('-r%d', cfg.figureDpi));
    
    figFileMatlab = fullfile(outDir, 'comparativo_pre_pos_bsgm.fig');
    savefig(fig, figFileMatlab);

    txtFile = fullfile(outDir, 'tabela_resultados_bsgm.txt');
    fid = fopen(txtFile, 'w');
    fprintf(fid, 'METRIC\tPRE_PRANDIAL\tPOST_PRANDIAL\tCOMBINED_VALUE\tUNIT_STATUS\n');
    fprintf(fid, 'Filter_Band_cpm\t%.2f-%.2f\t%.2f-%.2f\t-\tcpm\n', cfg.filterMinCpm, cfg.filterMaxCpm, cfg.filterMinCpm, cfg.filterMaxCpm);
    if cfg.hasBmi
        fprintf(fid, 'BMI_Adjustment\t%.1f\t%.1f\t-\tkg/m2\n', cfg.bmi, cfg.bmi);
    else
        fprintf(fid, 'BMI_Adjustment\tNot_Applied\tNot_Applied\t-\t-\n');
    end
    fprintf(fid, 'GA-RI\t%.6f\t%.6f\t-\tPre: %s | Post: %s\n', preResults.gari, postResults.gari, gari_status(preResults.gari), gari_status(postResults.gari));
    fprintf(fid, 'Mean_Amplitude\t%.6f\t%.6f\t-\tpT\n', preResults.meanAmplitude, postResults.meanAmplitude);
    fprintf(fid, 'Amplitude_P90\t%.6f\t%.6f\t-\tpT\n', postResults.amplitudeP90, preResults.meanAmplitude);
    fprintf(fid, 'Amplitude_Ratio\t-\t-\t%.6f\tratio\n', amplitudeRatio);
    fprintf(fid, 'ff-AR\t-\t-\t%.6f\tratio\n', fedFastedAmplitudeRatio);
    fclose(fid);
    fprintf('[OK] Results and figures saved to: %s\n\n', outDir);
end

%% ================= HELPER FUNCTIONS =================

function results = process_record(inputFile, cfg)
    M = read_numeric_compatible(inputFile);
    M = M(~all(isnan(M),2),:);
    M = M(:,~all(isnan(M),1));
    
    t = M(:,1);
    X = M(:,2:end) * cfg.picoTeslaFactor; % Retain the original values in pT
    [nSamples, nChannels] = size(X);
    
    timeIntervals = diff(t); timeIntervals = timeIntervals(isfinite(timeIntervals) & timeIntervals > 0);
    if isempty(timeIntervals), samplingRate = cfg.samplingRate; else, samplingRate = 1/median(timeIntervals); end
    
    centeredSignals = zeros(nSamples, nChannels);
    validChannels = true(1, nChannels);
    for ch = 1:nChannels
        x = X(:,ch);
        if mean(isfinite(x)) < 0.9 || max(x)-min(x) <= eps
            validChannels(ch) = false;
            continue;
        end
        centeredSignals(:,ch) = x - mean(x);
    end
    
    validChannelIndices = find(validChannels);
    globalSignal = mean(centeredSignals(:, validChannelIndices), 2);
    
    % =========================================================================
    % OPTION 1: Remove the linear trend (detrending)
    % Remove upward baseline drift/offset
    % =========================================================================
    detrendedGlobalSignal = detrend(globalSignal);

    windowDurationSeconds = 60;
    welchWindowSamples  = round(windowDurationSeconds * samplingRate);
    overlapSamples     = round((windowDurationSeconds / 2) * samplingRate);

    [powerSpectralDensity, frequencyHz] = pwelch(detrendedGlobalSignal, hann_compat(welchWindowSamples), overlapSamples, cfg.nfft, samplingRate);
    frequencyCpm = frequencyHz * 60;
    
    narrowBandMask = (frequencyCpm >= cfg.filterMinCpm & frequencyCpm <= cfg.filterMaxCpm);
    
    totalBandMinCpm = min(1.0, cfg.filterMinCpm);
    totalBandMaxCpm = max(6.0, cfg.filterMaxCpm);
    totalBandMask  = (frequencyCpm >= totalBandMinCpm & frequencyCpm <= totalBandMaxCpm);
    
    gari = sum(powerSpectralDensity(narrowBandMask)) / max(sum(powerSpectralDensity(totalBandMask)), eps);
    
    amplitudeWindowSamples = max(2, round(20 * samplingRate));
    windowCount = floor(nSamples / amplitudeWindowSamples);
    windowAmplitudes = zeros(windowCount, 1);
    for k = 1:windowCount
        segment = detrendedGlobalSignal((k-1)*amplitudeWindowSamples+1 : k*amplitudeWindowSamples);
        windowAmplitudes(k) = max(segment) - min(segment);
    end
    
    results.time = t;
    results.signal = detrendedGlobalSignal;
    results.frequencyCpm = frequencyCpm;
    results.powerSpectralDensity = powerSpectralDensity;
    results.gari = gari;
    results.meanAmplitude = mean(windowAmplitudes);
    results.amplitudeP90   = prctile(windowAmplitudes, 90);
end

function fig = plot_comparison(preResults, postResults, amplitudeRatio, fedFastedAmplitudeRatio, cfg)
    fig = figure('Color', 'w', 'Name', 'BSGM PRE vs POST comparison', 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
    
    globalMin = min([preResults.signal; postResults.signal]);
    globalMax = max([preResults.signal; postResults.signal]);
    margin = 0.05 * (globalMax - globalMin);
    if margin == 0, margin = 1; end
    yLimits = [globalMin - margin, globalMax + margin];

    subplot(2,2,1);
    plot(preResults.time / 60, preResults.signal, 'b', 'LineWidth', 1);
    title('Global pre-prandial signal');
    xlabel('Time (min)'); ylabel('Amplitude (pT)'); grid on;
    ylim(yLimits);
    
    subplot(2,2,2);
    plot(postResults.time / 60, postResults.signal, 'r', 'LineWidth', 1);
    title('Global post-prandial signal');
    xlabel('Time (min)'); ylabel('Amplitude (pT)'); grid on;
    ylim(yLimits);
    
    % Subplot 3: Power spectral density (PSD)
    subplot(2,2,3);
    hold on; grid on;
    
    plot(preResults.frequencyCpm, preResults.powerSpectralDensity, 'b', 'LineWidth', 1.5, 'DisplayName', 'Pre-prandial');
    plot(postResults.frequencyCpm, postResults.powerSpectralDensity, 'r', 'LineWidth', 1.5, 'DisplayName', 'Post-prandial');
    
    % Restrict peak selection to the normogastric range (2.0 - 4.0 cpm)
    normogastricMask = (preResults.frequencyCpm >= 2.0 & preResults.frequencyCpm <= 4.0);
    normogastricFrequencies = preResults.frequencyCpm(normogastricMask);
    
    if ~isempty(normogastricFrequencies)
        [prePeakPower, prePeakIndex] = max(preResults.powerSpectralDensity(normogastricMask));
        prePeakFrequency = normogastricFrequencies(prePeakIndex);
        
        [postPeakPower, postPeakIndex] = max(postResults.powerSpectralDensity(normogastricMask));
        postPeakFrequency = normogastricFrequencies(postPeakIndex);
        
        plot(prePeakFrequency, prePeakPower, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(postPeakFrequency, postPeakPower, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
        
        text(prePeakFrequency, prePeakPower * 1.04, sprintf('  Pre peak: %.2f cpm', prePeakFrequency), 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 8.5);
        text(postPeakFrequency, postPeakPower * 1.04, sprintf('  Post peak: %.2f cpm', postPeakFrequency), 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 8.5);
    end
        
    xlim([min(1.0, cfg.filterMinCpm - 0.2), max(6.0, cfg.filterMaxCpm + 0.5)]);
    
    filterBandMask = (preResults.frequencyCpm >= cfg.filterMinCpm & preResults.frequencyCpm <= cfg.filterMaxCpm);
    maxLocalPower = max([preResults.powerSpectralDensity(filterBandMask); postResults.powerSpectralDensity(filterBandMask)]);
    if maxLocalPower > 0, ylim([0, maxLocalPower * 1.25]); end
    
    title(sprintf('PSD (filter: %.1f - %.1f cpm | peak in normogastria)', cfg.filterMinCpm, cfg.filterMaxCpm));
    xlabel('Frequency (cpm)'); ylabel('PSD (pT^2/Hz)');
    legend('Location', 'northeast');
    
    subplot(2,2,4);
    axis off;
    
    preStatusText = gari_status(preResults.gari);
    postStatusText = gari_status(postResults.gari);
    
    if cfg.hasBmi
        bmiStatusText = sprintf('BMI adjustment: YES (%.1f kg/m^2)', cfg.bmi);
    else
        bmiStatusText = 'BMI adjustment: NO (raw values)';
    end
    
    summaryText = {
        'BSGM and EGG metrics:';
        sprintf('MGG filter band: %.2f - %.2f cpm', cfg.filterMinCpm, cfg.filterMaxCpm);
        bmiStatusText;
        '';
        sprintf('GA-RI Pre: %.3f %s', preResults.gari, preStatusText);
        sprintf('GA-RI Post: %.3f %s', postResults.gari, postStatusText);
        '';
        sprintf('Mean amplitude Pre: %.2f pT', preResults.meanAmplitude);
        sprintf('Mean amplitude Post: %.2f pT', postResults.meanAmplitude);
        sprintf('P90 amplitude Post: %.2f pT', postResults.amplitudeP90);
        '';
        sprintf('Amplitude Ratio: %.2f', amplitudeRatio);
        sprintf('ff-AR (Fed:Fasted P90): %.2f', fedFastedAmplitudeRatio)
    };
    
    text(0.05, 0.5, summaryText, 'Interpreter', 'none', 'FontSize', 10);
end

function statusText = gari_status(value)
    if value >= 0.25
        statusText = '(Normal >= 0.25)';
    else
        statusText = '(Abnormal < 0.25)';
    end
end

function M = read_numeric_compatible(filename)
    if exist('readmatrix','file') == 2
        try, M = readmatrix(filename); if isnumeric(M) && ~isempty(M), return; end; catch; end
    end
    try, M = dlmread(filename); return; catch; end
    imp = importdata(filename);
    if isstruct(imp) && isfield(imp,'data'), M = imp.data; else, M = imp; end
end

function w = hann_compat(n)
    if exist('hann','file') == 2, w = hann(n, 'periodic'); else, w = hanning(n); end
end