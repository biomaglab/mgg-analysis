function b_global_analysis(inputFile)
% B_GLOBAL_ANALYSIS
% Global processing and analysis of gastric bioelectrical/biomagnetic signals.
%
% Features:
% - Direct reading at the original picoTesla (pT) scale.
% - Interactive initial-time selection using a cursor/crosshair (ginput).
% - Preliminary per-channel quality control (data, plateau, constant signal, RMS).
% - Global spectral analysis using FFT, Welch, and RSA (Running Spectral Analysis).
% - Instability Coefficient (IC) calculation based on the GastroPy model.
% - Automatic data export to .DAT, .TXT, .FIG, and .PNG files.

    if nargin < 1 || isempty(inputFile)
        [file, folder] = uigetfile({'*.dat;*.txt;*.csv;*.tsv','Numeric files'; '*.*','All files'}, ...
            'Select the post-pipeline file');
        if isequal(file,0)
            fprintf('No file selected.\n');
            return;
        end
        inputFile = fullfile(folder,file);
    end

    %% Pipeline configuration
    cfg.fsPadrao = 20.0;
    cfg.freqMinHz = 0.015;
    cfg.freqMaxHz = 0.200;
    cfg.tInicioSeg = 50.0; % Default initial time to discard
    cfg.janelaWelchSeg = 300.0;
    cfg.sobreposicaoWelch = 0.75;
    cfg.nfftWelch = 32768;
    cfg.nfftFFT = 32768;
    cfg.janelaRSASeg = 120.0;
    cfg.passoRSASeg = 15.0;
    cfg.nfftRSA = 16384;
    cfg.janelaAmplitudeSeg = 20.0;
    cfg.fracaoMinValida = 0.90;
    cfg.duracaoPlatoMaxSeg = 2.0;
    cfg.tolPlatoRelativa = 1e-10;
    cfg.fatorRMSQuaseNulo = 0.01;
    cfg.dpiFiguras = 180;
    cfg.fatorPicoTesla = 1.0; % Unit factor: the signal is already in pT

    BRADI = [1 2]; NORMO = [2 4]; TAQUI = [4 9]; INTESTINAL = [9 12];
    bandNames = {'Bradygastria','Normogastria','Tachygastria','Intestinal band'};

    %% Read data
    M = read_numeric_compatible(inputFile);
    M = M(~all(isnan(M),2),:);
    M = M(:,~all(isnan(M),1));
    if size(M,2) < 2 || size(M,1) < 2
        error('Invalid file: expected a time vector in column 1 and at least one channel.');
    end
    tempoOrig = M(:,1);
    XorigBruto = M(:,2:end) * cfg.fatorPicoTesla;

    dt = diff(tempoOrig);
    dt = dt(isfinite(dt) & dt > 0);
    if isempty(dt)
        Fs = cfg.fsPadrao;
    else
        Fs = 1/median(dt);
        if abs(Fs-cfg.fsPadrao) <= 0.01*cfg.fsPadrao
            Fs = cfg.fsPadrao;
        end
    end
    if ~isfinite(Fs) || Fs <= 0, error('Invalid sampling frequency.'); end

    fprintf('\nFile: %s\n', inputFile);
    fprintf('Detected Fs: %.8g Hz\n', Fs);
    fprintf('Original total duration: %.3f s (%.3f min)\n', tempoOrig(end)-tempoOrig(1), (tempoOrig(end)-tempoOrig(1))/60);

    %% 1. Preview and time cropping using the cursor (crosshair)
    [tClicado, figPrevia] = plot_raw_signal_preview_with_crosshair(tempoOrig, XorigBruto, cfg.tInicioSeg);
    
    if isfinite(tClicado)
        cfg.tInicioSeg = tClicado;
    end

    fprintf('\n-----------------------------------------------------------------\n');
    fprintf('Time selected visually: %.2f s\n', cfg.tInicioSeg);
    fprintf('-----------------------------------------------------------------\n');

    r = input(sprintf('Press Enter to confirm %.2f s or enter a new value (0 = no crop): ', cfg.tInicioSeg), 's');
    if ~isempty(strtrim(r))
        valInput = str2double(r);
        if isfinite(valInput) && valInput >= 0
            cfg.tInicioSeg = valInput;
        end
    end

    if isgraphics(figPrevia) && isvalid(figPrevia)
        close(figPrevia);
    end

    % Apply the initial time crop
    tempo = tempoOrig;
    Xorig = XorigBruto;
    if cfg.tInicioSeg > 0
        idxInicio = find(tempoOrig >= tempoOrig(1) + cfg.tInicioSeg, 1, 'first');
        if ~isempty(idxInicio) && idxInicio < numel(tempoOrig)
            tempo = tempoOrig(idxInicio:end);
            Xorig = XorigBruto(idxInicio:end, :);
            fprintf('--> Signal cropped. Analyzing from t = %.2f s (%.2f min) to t = %.2f s.\n', tempo(1), tempo(1)/60, tempo(end));
        else
            warning('The specified crop time exceeds the signal duration. Cropping was skipped.');
        end
    else
        fprintf('--> No crop applied. Analyzing the full signal starting at t = %.2f s.\n', tempo(1));
    end

    [nAmostras, nCanais] = size(Xorig);
    fprintf('Usable samples: %d | Channels: %d\n', nAmostras, nCanais);
    fprintf('Usable duration: %.3f s (%.3f min)\n', tempo(end)-tempo(1), (tempo(end)-tempo(1))/60);
    fprintf('Amplitude scale: Retained at input scale (picoTesla - pT)\n');

    %% Set the analysis frequency band
    fprintf('\nThe routine does NOT reapply the band-pass filter.\n');
    fprintf('Enter only the spectral range to analyze.\n');
    r = input(sprintf('Lower frequency [%.6g Hz = %.3f cpm]: ', cfg.freqMinHz, cfg.freqMinHz*60), 's');
    if ~isempty(strtrim(r)), cfg.freqMinHz = str2double(r); end
    r = input(sprintf('Upper frequency [%.6g Hz = %.3f cpm]: ', cfg.freqMaxHz, cfg.freqMaxHz*60), 's');
    if ~isempty(strtrim(r)), cfg.freqMaxHz = str2double(r); end
    if ~(isfinite(cfg.freqMinHz)&&isfinite(cfg.freqMaxHz)&&0<cfg.freqMinHz&&cfg.freqMinHz<cfg.freqMaxHz&&cfg.freqMaxHz<Fs/2)
        error('Invalid spectral band.');
    end
    cfg.cpmMin = cfg.freqMinHz*60;
    cfg.cpmMax = cfg.freqMaxHz*60;
    fprintf('Analyzed band: %.4f-%.4f cpm (%.6g-%.6g Hz).\n', cfg.cpmMin, cfg.cpmMax, cfg.freqMinHz, cfg.freqMaxHz);

    %% Per-channel technical quality control
    Xtrat = nan(size(Xorig));
    fracaoValida = nan(1,nCanais); rmsCanal = nan(1,nCanais); maiorPlatoSeg = nan(1,nCanais);
    falhaDados = false(1,nCanais); falhaConstante = false(1,nCanais); falhaPlato = false(1,nCanais); falhaRMS = false(1,nCanais);

    for ch=1:nCanais
        x = Xorig(:,ch);
        validos = isfinite(x);
        fracaoValida(ch) = mean(validos);
        if nnz(validos)<2 || fracaoValida(ch)<cfg.fracaoMinValida
            falhaDados(ch)=true; continue;
        end
        if ~all(validos)
            idx=(1:nAmostras)';
            x(~validos)=interp1(idx(validos),x(validos),idx(~validos),'linear','extrap');
        end
        x=x-mean(x);
        ampTotal=max(x)-min(x);
        if ~isfinite(ampTotal)||ampTotal<=eps
            falhaConstante(ch)=true; continue;
        end
        tolerancia=max(eps(ampTotal),ampTotal*cfg.tolPlatoRelativa);
        repetido=abs(diff(x))<=tolerancia;
        maiorSeq=longest_true_run(repetido);
        maiorPlatoSeg(ch)=(maiorSeq+1)/Fs;
        falhaPlato(ch)=maiorPlatoSeg(ch)>=cfg.duracaoPlatoMaxSeg;
        Xtrat(:,ch)=x;
        rmsCanal(ch)=sqrt(mean(x.^2));
    end

    candidatos=find(~falhaDados & ~falhaConstante & ~falhaPlato & isfinite(rmsCanal));
    if isempty(candidatos), error('No technically valid channels were found.'); end
    medRMS=median(rmsCanal(candidatos));
    limiarRMS=cfg.fatorRMSQuaseNulo*medRMS;
    falhaRMS(candidatos)=rmsCanal(candidatos)<limiarRMS;
    alertaRMS=falhaRMS;
    falhaAuto=falhaDados|falhaConstante|falhaPlato;

    fprintf('\nPreliminary technical quality control:\n');
    for ch=1:nCanais
        avisos={};
        if falhaDados(ch), avisos{end+1}='insufficient data'; end
        if falhaConstante(ch), avisos{end+1}='constant signal'; end
        if falhaPlato(ch), avisos{end+1}='prolonged plateau'; end
        if alertaRMS(ch), avisos{end+1}='near-zero RMS - WARNING'; end
        if isempty(avisos), txt='no technical warnings'; else txt=strjoin_compatible(avisos,'; '); end
        fprintf('  Channel %d: %s\n',ch,txt);
    end

    %% Visual channel inspection
    figInspecao=plot_channel_inspection(tempo,Xorig,falhaDados,falhaConstante,falhaPlato,alertaRMS);
    fprintf('\nInspect the channel quality figure. Close the window to continue.\n');
    waitfor(figInspecao);

    fprintf('\nChannels with objective failures (data/constant/plateau) will be excluded automatically.\n');
    fprintf('Low RMS is a warning only and does NOT cause automatic exclusion.\n');
    resposta=input('Additional channels to exclude [e.g. 2 5 7 | Enter = none]: ','s');
    exclusaoManual=parse_channel_list(resposta,nCanais);

    canaisValidosIdx=find(~falhaAuto & ~ismember(1:nCanais,exclusaoManual));
    canaisExcluidosIdx=setdiff(1:nCanais,canaisValidosIdx);
    if isempty(canaisValidosIdx), error('No valid channels remain for global analysis.'); end
    fprintf('\nChannels used in the global average: %s\n',num2str(canaisValidosIdx));
    if isempty(canaisExcluidosIdx), fprintf('Excluded channels: none\n'); else fprintf('Excluded channels: %s\n',num2str(canaisExcluidosIdx)); end

    %% Global mean and median signals
    sinalGlobal=mean(Xtrat(:,canaisValidosIdx),2);
    sinalGlobal=sinalGlobal-mean(sinalGlobal);
    sinalGlobalMediana=median(Xtrat(:,canaisValidosIdx),2);
    sinalGlobalMediana=sinalGlobalMediana-mean(sinalGlobalMediana);

    %% Global FFT and Welch
    fftInfo=calc_fft(sinalGlobal,Fs,cfg.cpmMin,cfg.cpmMax,cfg.nfftFFT);
    welchInfo=calc_welch(sinalGlobal,Fs,cfg);
    classeGlobal=classify_frequency(welchInfo.picoCpm);

    potenciaTotal=trapz(welchInfo.fLocalCpm,welchInfo.pLocal);
    potFaixas=[integrate_band(welchInfo.fCpm,welchInfo.pxx,BRADI), ...
               integrate_band(welchInfo.fCpm,welchInfo.pxx,NORMO), ...
               integrate_band(welchInfo.fCpm,welchInfo.pxx,TAQUI), ...
               integrate_band(welchInfo.fCpm,welchInfo.pxx,INTESTINAL)];
    if sum(potFaixas)>0, fracaoPotencia=100*potFaixas/sum(potFaixas); else, fracaoPotencia=nan(1,4); end
    fwhm=calc_fwhm(welchInfo.fLocalCpm,welchInfo.pLocal,welchInfo.idxPicoLocal);
    entropia=spectral_entropy(welchInfo.pLocal);

    %% Calculate the Instability Coefficient (IC)
    icGlobalInfo = calc_instability_coefficient(sinalGlobal, Fs);

    janelaAmp=max(2,round(cfg.janelaAmplitudeSeg*Fs));
    nJan=floor(nAmostras/janelaAmp); ampPPJan=nan(nJan,1);
    for k=1:nJan
        seg=sinalGlobal((k-1)*janelaAmp+1:k*janelaAmp);
        ampPPJan(k)=max(seg)-min(seg);
    end
    if isempty(ampPPJan), ampPPMedia=NaN; ampPPMediana=NaN; else, ampPPMedia=mean(ampPPJan); ampPPMediana=median(ampPPJan); end
    ampPPGlobal=max(sinalGlobal)-min(sinalGlobal);
    rmsGlobal=sqrt(mean(sinalGlobal.^2));

    %% RSA (Running Spectral Analysis)
    rsa=calc_rsa(sinalGlobal,tempo,Fs,cfg);

    %% Individual per-channel analysis
    controle=struct([]);
    limiarConcordancia=0.25;
    for ch=1:nCanais
        controle(ch).canal=ch;
        controle(ch).valido_pct=100*fracaoValida(ch);
        controle(ch).rms=rmsCanal(ch);
        if ismember(ch,canaisValidosIdx)
            fi=calc_fft(Xtrat(:,ch),Fs,cfg.cpmMin,cfg.cpmMax,cfg.nfftFFT);
            wi=calc_welch(Xtrat(:,ch),Fs,cfg);
            icChInfo = calc_instability_coefficient(Xtrat(:,ch), Fs);
            
            controle(ch).fft_cpm=fi.picoCpm;
            controle(ch).welch_cpm=wi.picoCpm;
            controle(ch).classe_fft=classify_frequency(fi.picoCpm);
            controle(ch).classe_welch=classify_frequency(wi.picoCpm);
            controle(ch).delta=abs(fi.picoCpm-wi.picoCpm);
            if controle(ch).delta<=limiarConcordancia
                controle(ch).concordancia='AGREE';
            elseif strcmp(controle(ch).classe_fft,controle(ch).classe_welch)
                controle(ch).concordancia='SAME_BAND_DIFFERENT_PEAKS';
            else
                controle(ch).concordancia='DIFFERENT_BANDS';
            end
            controle(ch).entropia=spectral_entropy(wi.pLocal);
            controle(ch).ic=icChInfo.IC;
            if alertaRMS(ch), controle(ch).status='VALID_WITH_RMS_WARNING'; else, controle(ch).status='VALID'; end
        else
            controle(ch).fft_cpm=NaN; controle(ch).welch_cpm=NaN;
            controle(ch).classe_fft='NOT_ANALYZED'; controle(ch).classe_welch='NOT_ANALYZED';
            controle(ch).delta=NaN; controle(ch).concordancia='NOT_ANALYZED'; controle(ch).entropia=NaN;
            controle(ch).ic=NaN;
            if ismember(ch,exclusaoManual), controle(ch).status='MANUAL_EXCLUSION_CONFIRMED_FAILURE';
            elseif falhaDados(ch), controle(ch).status='DATA_FAILURE';
            elseif falhaConstante(ch), controle(ch).status='CONSTANT_SIGNAL';
            elseif falhaPlato(ch), controle(ch).status='PROLONGED_PLATEAU';
            else, controle(ch).status='EXCLUDED'; end
        end
    end

    %% Output directory
    [folder, base, ~] = fileparts(inputFile);
    cleanBase = strrep(base, '_original_dec20Hz', '');
    cleanBase = strrep(cleanBase, '_pipeline_sem_ICA', '');
    cleanBase = strrep(cleanBase, '_analise_global', '');

    if ~isempty(strfind(lower(cleanBase), 'pos'))
        suffix = '(pos)';
    else
        suffix = '(pre)';
    end

    outDir = fullfile(folder, [cleanBase suffix]);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    %% Generate figures
    figs(1) = plot_stacked_channels(tempo, Xtrat, canaisValidosIdx);
    figs(2) = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
    plot(tempo/60, sinalGlobal, 'LineWidth', 0.9); grid on; 
    xlabel('Time (min)'); ylabel('Amplitude (pT)'); 
    title(sprintf('Global mean trace - %d valid channels (Starting at %.1f s)', numel(canaisValidosIdx), tempo(1)));

    figs(3) = plot_fft_welch(fftInfo, welchInfo, cfg, classeGlobal, fwhm);
    figs(4) = plot_rsa(tempo, sinalGlobal, rsa, cfg);
    figs(5) = plot_classification(fracaoPotencia, rsa.prevalencia);

    figureBaseNames = {'01_sinais_bandpass', ...
                 '02_tracado_global', ...
                 '03_fft_welch_global', ...
                 '04_rsa_global_e_tracado', ...
                 '05_classificacao_fisiologica'};

    for k = 1:numel(figs)
        if isgraphics(figs(k)) && isvalid(figs(k))
            print(figs(k), fullfile(outDir, [figureBaseNames{k} '.png']), '-dpng', sprintf('-r%d', cfg.dpiFiguras));
            savefig(figs(k), fullfile(outDir, [figureBaseNames{k} '.fig']));
        end
    end

    %% Export numeric data to .DAT files
    % A) Global time-domain trace in .DAT
    matrizSinaisDAT = [tempo, tempo/60, sinalGlobal, sinalGlobalMediana];
    fileSinaisDAT = fullfile(outDir, 'sinal_global_tempo.dat');
    fidDAT = fopen(fileSinaisDAT, 'w');
    fprintf(fidDAT, 'Time_s\tTime_min\tSignal_Mean_pT\tSignal_Median_pT\n');
    fclose(fidDAT);
    dlmwrite(fileSinaisDAT, matrizSinaisDAT, '-append', 'delimiter', '\t', 'precision', '%.8f');

    % B) Welch and FFT frequency spectrum in .DAT
    matrizEspectroDAT = [welchInfo.fLocalCpm/60, welchInfo.fLocalCpm, welchInfo.pLocal, fftInfo.aLocal];
    fileEspectroDAT = fullfile(outDir, 'espectro_welch_fft.dat');
    fidDAT = fopen(fileEspectroDAT, 'w');
    fprintf(fidDAT, 'Freq_Hz\tFreq_cpm\tPSD_Welch_pT2_Hz\tFFT_Amplitude_pT\n');
    fclose(fidDAT);
    dlmwrite(fileEspectroDAT, matrizEspectroDAT, '-append', 'delimiter', '\t', 'precision', '%.8f');

    % C) Channel quality-control table in .DAT
    fileControleDAT = fullfile(outDir, 'controle_canais.dat');
    fidDAT = fopen(fileControleDAT, 'w');
    fprintf(fidDAT, 'Channel\tValid_pct\tRMS_pT\tFFT_cpm\tFFT_Class\tWelch_cpm\tWelch_Class\tDelta\tAgreement\tEntropy\tIC\tStatus\n');
    for ch=1:nCanais
        fprintf(fidDAT, '%d\t%.6f\t%.10g\t%.8f\t%s\t%.8f\t%s\t%.8f\t%s\t%.8f\t%.8f\t%s\n', ...
            controle(ch).canal, controle(ch).valido_pct, controle(ch).rms, controle(ch).fft_cpm, controle(ch).classe_fft, ...
            controle(ch).welch_cpm, controle(ch).classe_welch, controle(ch).delta, controle(ch).concordancia, controle(ch).entropia, controle(ch).ic, controle(ch).status);
    end
    fclose(fidDAT);

    % D) Filtered channel signals in .DAT
    matrizCanaisDAT = [tempo, Xtrat];
    fileCanaisDAT = fullfile(outDir, 'sinais_canais_filtrados.dat');
    fidDAT = fopen(fileCanaisDAT, 'w');
    headersCanais = sprintf('Time_s\t');
    for ch = 1:nCanais
        headersCanais = [headersCanais, sprintf('Channel_%d_pT\t', ch)]; %#ok<AGROW>
    end
    headersCanais = [strtrim(headersCanais), '\n'];
    fprintf(fidDAT, headersCanais);
    fclose(fidDAT);
    dlmwrite(fileCanaisDAT, matrizCanaisDAT, '-append', 'delimiter', '\t', 'precision', '%.8f');

    %% TXT report
    fid=fopen(fullfile(outDir,'relatorio_analise_global.txt'),'w');
    fprintf(fid,'POST-PIPELINE GLOBAL GASTRIC ANALYSIS\n===================================\n\n');
    fprintf(fid,'File: %s\nFs: %.10g Hz\nSamples: %d\nUsable duration: %.6f s (Initial time discarded: %.1f s)\n',inputFile,Fs,nAmostras,tempo(end)-tempo(1),cfg.tInicioSeg);
    fprintf(fid,'Amplitude scale: picoTesla (pT) [direct/no conversion]\n');
    fprintf(fid,'Analyzed spectral band: %.6g-%.6g cpm (%.8g-%.8g Hz)\n',cfg.cpmMin,cfg.cpmMax,cfg.freqMinHz,cfg.freqMaxHz);
    fprintf(fid,'Valid channels: %s\nExcluded channels: %s\n\n',num2str(canaisValidosIdx),num2str(canaisExcluidosIdx));
    fprintf(fid,'MAIN GLOBAL RESULT\n');
    fprintf(fid,'Global FFT peak: %.8f cpm\n',fftInfo.picoCpm);
    fprintf(fid,'Welch dominant frequency: %.8f cpm (%.8f Hz)\n',welchInfo.picoCpm,welchInfo.picoCpm/60);
    fprintf(fid,'Global classification: %s\n',classeGlobal);
    fprintf(fid,'Global RMS (pT): %.10g\nGlobal peak-to-peak amplitude (pT): %.10g\n',rmsGlobal,ampPPGlobal);
    fprintf(fid,'Mean window peak-to-peak amplitude over %.1f s (pT): %.10g\nMedian window peak-to-peak amplitude (pT): %.10g\n',cfg.janelaAmplitudeSeg,ampPPMedia,ampPPMediana);
    fprintf(fid,'Total power %.3f-%.3f cpm (pT^2): %.10g\nPower at Welch peak (pT^2/Hz): %.10g\n',cfg.cpmMin,cfg.cpmMax,potenciaTotal,welchInfo.potenciaPico);
    fprintf(fid,'FWHM: %.8f cpm | bounds %.8f to %.8f cpm\n',fwhm(1),fwhm(2),fwhm(3));
    fprintf(fid,'Normalized spectral entropy: %.8f\n',entropia);
    fprintf(fid,'Instability Coefficient (IC): %.8f (Cycles: %d, Mean duration: %.2f s, Std. dev.: %.2f s)\n\n', ...
        icGlobalInfo.IC, icGlobalInfo.nCiclos, icGlobalInfo.meanDur, icGlobalInfo.stdDur);
    
    fprintf(fid,'POWER FRACTION BY BAND (%%)\n');
    for k=1:4, fprintf(fid,'%s: %.6f\n',bandNames{k},fracaoPotencia(k)); end
    fprintf(fid,'\nRSA TEMPORAL PREVALENCE (%%)\n');
    for k=1:4, fprintf(fid,'%s: %.6f\n',bandNames{k},rsa.prevalencia(k)); end
    fprintf(fid,'\nINTERNAL PER-CHANNEL QUALITY CONTROL\nChannel\tValid_pct\tRMS(pT)\tFFT_cpm\tFFT_Class\tWelch_cpm\tWelch_Class\tDelta\tAgreement\tEntropy\tIC\tStatus\n');
    for ch=1:nCanais
        fprintf(fid,'%d\t%.6f\t%.10g\t%.8f\t%s\t%.8f\t%s\t%.8f\t%s\t%.8f\t%.8f\t%s\n', ...
            controle(ch).canal,controle(ch).valido_pct,controle(ch).rms,controle(ch).fft_cpm,controle(ch).classe_fft, ...
            controle(ch).welch_cpm,controle(ch).classe_welch,controle(ch).delta,controle(ch).concordancia,controle(ch).entropia,controle(ch).ic,controle(ch).status);
    end
    fclose(fid);

    save(fullfile(outDir,'analise_global_completa.mat'),'tempo','Fs','canaisValidosIdx','canaisExcluidosIdx','sinalGlobal','sinalGlobalMediana', ...
        'fftInfo','welchInfo','classeGlobal','potenciaTotal','potFaixas','fracaoPotencia','fwhm','entropia','icGlobalInfo','ampPPGlobal','ampPPMedia','ampPPMediana','rmsGlobal','rsa','controle','cfg');

    fprintf('\n=============================================\nGLOBAL ANALYSIS COMPLETE\n=============================================\n');
    fprintf('Valid channels: %s\n',num2str(canaisValidosIdx));
    fprintf('Welch dominant frequency: %.4f cpm\n',welchInfo.picoCpm);
    fprintf('Classification: %s\n',classeGlobal);
    fprintf('Results and .DAT files saved to:\n  %s\n',outDir);

    %% Optional individual plots
    while true
        r=input('\nPlot individual channel with signal + FFT + Welch? [number | 0=exit]: ','s');
        if isempty(r) || strcmp(strtrim(r),'0'), break; end
        ch=str2double(r);
        if ~isfinite(ch)||ch<1||ch>nCanais||~ismember(ch,canaisValidosIdx)
            fprintf('Invalid or excluded channel.\n'); continue;
        end
        plot_individual_channel(ch,tempo,Xtrat(:,ch),Fs,cfg);
    end
end

%% ================= COMPATIBLE DATA-READING HELPER =================
function M = read_numeric_compatible(filename)
    if exist('readmatrix', 'file') == 2
        try
            M = readmatrix(filename);
            if isnumeric(M) && ~isempty(M)
                return;
            end
        catch
        end
    end

    try
        M = dlmread(filename);
        return;
    catch
    end

    imp = importdata(filename);
    if isstruct(imp) && isfield(imp, 'data')
        M = imp.data;
    else
        M = imp;
    end

    if ~isnumeric(M) || isempty(M)
        error('Could not read a valid numeric matrix from the provided file.');
    end
end

%% ================= CROSSHAIR HELPER FOR TIME-CROP SELECTION =================
function [tCortado, fig] = plot_raw_signal_preview_with_crosshair(t, X, tSugerido)
    nCanais = size(X, 2);
    fig = figure('Color', 'w', 'Name', 'CROSSHAIR TIME-CROP SELECTION', 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9]);
    hold on;
    
    amps = zeros(1, nCanais);
    for k = 1:nCanais
        x = X(:,k);
        amps(k) = percentile_compat(x(~isnan(x)), 95) - percentile_compat(x(~isnan(x)), 5);
    end
    esc = median(amps(isfinite(amps) & amps > 0));
    if isempty(esc) || ~isfinite(esc) || esc <= 0, esc = 1; end
    off = 1.6 * esc;
    
    bases = zeros(1, nCanais);
    for ch = 1:nCanais
        x = X(:,ch) - mean(X(:,ch), 'omitnan');
        base = (nCanais - ch) * off;
        bases(ch) = base;
        plot(t, x + base, 'LineWidth', 0.75);
    end
    
    yl = ylim;
    line([tSugerido tSugerido], yl, 'Color', [0.7 0.7 0.7], 'LineStyle', ':', 'LineWidth', 1.5);
    
    [basesOrd, ordemTicks] = sort(bases, 'ascend');
    rot = arrayfun(@(c) sprintf('Ch %d', c), (1:nCanais), 'UniformOutput', false);
    
    set(gca, 'YTick', basesOrd, 'YTickLabel', rot(ordemTicks));
    xlabel('Time (seconds)');
    ylabel('Channels (pT)');
    title({'PREVIEW - All channels (complete signal)', ...
           'CLICK WITH THE CROSSHAIR on the plot to set the start time (or press ENTER to use 0 s)'});
    grid on;
    xlim([t(1) t(end)]);
    drawnow;

    try
        [xClick, ~, button] = ginput(1);
        if ~isempty(xClick) && ~isempty(button) && button == 1
            tCortado = max(0, xClick);
            line([tCortado tCortado], yl, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            text(tCortado + 2, yl(2)*0.9, sprintf('Selected crop: %.2f s', tCortado), ...
                'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10, 'BackgroundColor', 'w');
            drawnow;
            pause(0.5);
        else
            tCortado = 0.0;
        end
    catch
        tCortado = NaN;
    end
end

%% ================= IC HELPER FUNCTIONS =================
function infoIC = calc_instability_coefficient(x, Fs)
    x = x(:);
    if numel(x) < 2 || all(x == 0)
        infoIC = struct('IC', NaN, 'nCiclos', 0, 'meanDur', NaN, 'stdDur', NaN);
        return;
    end
    
    fase = angle(hilbert(x));
    transicoes = find(diff(fase) < -1.8*pi);
    
    if numel(transicoes) < 2
        infoIC = struct('IC', NaN, 'nCiclos', numel(transicoes), 'meanDur', NaN, 'stdDur', NaN);
        return;
    end
    
    duracoesSeg = diff(transicoes) / Fs;
    mediaDur = mean(duracoesSeg);
    stdDur = std(duracoesSeg);
    
    if mediaDur > 0
        ic = stdDur / mediaDur;
    else
        ic = NaN;
    end
    
    infoIC = struct('IC', ic, 'nCiclos', numel(duracoesSeg), 'meanDur', mediaDur, 'stdDur', stdDur);
end

%% ================= OTHER HELPER FUNCTIONS =================
function n=longest_true_run(v)
    n=0; atual=0;
    for i=1:numel(v)
        if v(i), atual=atual+1; n=max(n,atual); else, atual=0; end
    end
end

function lista=parse_channel_list(txt,nMax)
    txt=strrep(txt,',',' '); txt=strrep(txt,'[',''); txt=strrep(txt,']','');
    if isempty(strtrim(txt)), lista=[]; return; end
    lista=unique(sscanf(txt,'%d').');
    if any(lista<1|lista>nMax), error('Channel indices are outside the allowed range.'); end
end

function fig=plot_channel_inspection(t,X,fd,fc,fp,fr)
    n=size(X,2); fig=figure('Color','w','Name','GLOBAL QUALITY INSPECTION','NumberTitle','off','Units','normalized','OuterPosition',[0 0 1 1]); hold on;
    bases=zeros(1,n); rot=cell(1,n);
    for ch=1:n
        x=X(:,ch); good=isfinite(x);
        if nnz(good)>=2
            if ~all(good), idx=(1:numel(x))'; x(~good)=interp1(idx(good),x(good),idx(~good),'linear','extrap'); end
            x=x-mean(x); esc=percentile_compat(x,95)-percentile_compat(x,5); if ~isfinite(esc)||esc<=eps, esc=std(x); end
            if isfinite(esc)&&esc>eps, x=x/esc; else, x(:)=NaN; end
        else, x(:)=NaN; end
        base=(n-ch)*2; bases(ch)=base; plot(t/60,x+base,'LineWidth',0.75);
        av={}; if fd(ch),av{end+1}='DATA';end; if fc(ch),av{end+1}='CONSTANT';end; if fp(ch),av{end+1}='PLATEAU';end; if fr(ch),av{end+1}='LOW RMS';end
        if isempty(av), rot{ch}=sprintf('Ch %d',ch); else, rot{ch}=sprintf('Ch %d [%s]',ch,strjoin_compatible(av,' | ')); end
    end
    [basesOrd, ordemTicks] = sort(bases, 'ascend');
    rotOrd = rot(ordemTicks);
    set(gca,'YTick',basesOrd,'YTickLabel',rotOrd); xlabel('Time (min)'); ylabel('Channels (normalized for inspection only)');
    title({'CHANNEL QUALITY INSPECTION (POST-CROP)','Close this window to specify channels to exclude'}); grid on; xlim([t(1) t(end)]/60);
end

function fig=plot_stacked_channels(t,X,idx)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); hold on; n=numel(idx); amps=zeros(1,n);
    for k=1:n, x=X(:,idx(k)); amps(k)=percentile_compat(x,95)-percentile_compat(x,5); end
    esc=median(amps(isfinite(amps)&amps>0)); if isempty(esc), esc=max(std(X(:,idx),0,1)); end; if ~isfinite(esc)||esc<=0,esc=1;end
    off=1.6*esc; bases=zeros(1,n);
    for k=1:n, bases(k)=(n-k)*off; plot(t/60,X(:,idx(k))+bases(k),'LineWidth',0.75); end
    rotEmp = arrayfun(@(c)sprintf('Ch %d',c),idx,'UniformOutput',false);
    [basesOrd, ordemTicks] = sort(bases, 'ascend');
    rotEmpOrd = rotEmp(ordemTicks);
    set(gca,'YTick',basesOrd,'YTickLabel',rotEmpOrd); xlabel('Time (min)'); ylabel('Filtered channels (pT)'); title(sprintf('Stacked post-pipeline signals - %d channels used',n)); grid on; xlim([t(1) t(end)]/60);
end

function info=calc_fft(x,Fs,cpmMin,cpmMax,nfftMin)
    x=x(:); nfft=max(nfftMin,2^nextpow2(numel(x))); w=hann_compat(numel(x)); Y=fft(x.*w,nfft); npos=floor(nfft/2)+1; Y=Y(1:npos); fHz=(0:npos-1)'*Fs/nfft;
    amp=abs(Y)/max(sum(w),eps); if numel(amp)>2, amp(2:end-1)=2*amp(2:end-1); end
    fCpm=fHz*60; m=fCpm>=cpmMin & fCpm<=cpmMax; fl=fCpm(m); al=amp(m); [~,ix]=max(al); pico=refine_peak(fl,al,ix); ap=interp1(fl,al,pico,'linear');
    info=struct('fHz',fHz,'fCpm',fCpm,'amplitude',amp,'fLocalCpm',fl,'aLocal',al,'picoCpm',pico,'amplitudePico',ap);
end

function info=calc_welch(x,Fs,cfg)
    n=numel(x); wn=min(round(cfg.janelaWelchSeg*Fs),n); wn=max(min(wn,n),min(16,n)); w=hann_compat(wn); nov=min(wn-1,floor(cfg.sobreposicaoWelch*wn)); nfft=max(cfg.nfftWelch,2^nextpow2(wn));
    [pxx,fHz]=pwelch(x,w,nov,nfft,Fs,'psd'); fCpm=fHz*60; m=fCpm>=cfg.cpmMin & fCpm<=cfg.cpmMax; fl=fCpm(m); pl=pxx(m); [~,ix]=max(pl); pico=refine_peak(fl,pl,ix); pp=interp1(fl,pl,pico,'linear');
    info=struct('fHz',fHz,'fCpm',fCpm,'pxx',pxx,'fLocalCpm',fl,'pLocal',pl,'idxPicoLocal',ix,'picoCpm',pico,'potenciaPico',pp,'janelaN',wn,'noverlap',nov,'nfft',nfft);
end

function pico=refine_peak(f,y,ix)
    pico=f(ix); if ix>1 && ix<numel(y) && numel(f)>=2
        y1=log(max(y(ix-1),realmin)); y2=log(max(y(ix),realmin)); y3=log(max(y(ix+1),realmin)); den=y1-2*y2+y3;
        if isfinite(den)&&abs(den)>eps, d=0.5*(y1-y3)/den; d=max(-1,min(1,d)); pico=f(ix)+d*(f(2)-f(1)); end
    end
end

function val=integrate_band(f,p,band)
    m=f>=band(1)&f<=band(2); if nnz(m)>=2,val=trapz(f(m),p(m));elseif nnz(m)==1,val=p(m);else,val=0;end
end

function f=calc_fwhm(freq,p,ix)
    meia=p(ix)/2; fe=NaN; fd=NaN;
    iesq=find(p(1:ix)<=meia,1,'last'); if ~isempty(iesq)&&iesq<ix, fe=interp_cross(freq(iesq),p(iesq),freq(iesq+1),p(iesq+1),meia); end
    idr0=find(p(ix:end)<=meia,1,'first'); if ~isempty(idr0), idr=ix+idr0-1; if idr>ix, fd=interp_cross(freq(idr-1),p(idr-1),freq(idr),p(idr),meia); end, end
    if isfinite(fe)&&isfinite(fd)&&fd>fe, largura=fd-fe; else, largura=NaN; end
    f=[largura fe fd];
end

function x=interp_cross(f1,p1,f2,p2,a), if abs(p2-p1)>eps,x=f1+(a-p1)*(f2-f1)/(p2-p1);else,x=(f1+f2)/2;end, end

function h=spectral_entropy(p)
    p=p(isfinite(p)&p>=0); if numel(p)<2||sum(p)<=0,h=NaN;return;end; q=p/sum(p);q=q(q>0);h=-sum(q.*log(q))/log(numel(p));h=max(0,min(1,h));
end

function rsa=calc_rsa(x,t,Fs,cfg)
    n=numel(x); wn=min(round(cfg.janelaRSASeg*Fs),n); wn=max(min(wn,n),min(16,n)); passo=max(1,round(cfg.passoRSASeg*Fs));
    if wn>=n, inicios=1; else, inicios=1:passo:(n-wn+1); end
    w=hann_compat(wn); nov=floor(0.5*wn); nfft=max(cfg.nfftRSA,2^nextpow2(wn));
    potNorm=[]; freq=[]; tm=nan(numel(inicios),1); fd=nan(numel(inicios),1); cl=nan(numel(inicios),1);
    for k=1:numel(inicios)
        ii=inicios(k):(inicios(k)+wn-1); [pxx,fHz]=pwelch(x(ii),w,nov,nfft,Fs,'psd'); fc=fHz*60; m=fc>=cfg.cpmMin&fc<=cfg.cpmMax; fl=fc(m); pl=pxx(m);
        if isempty(potNorm), potNorm=nan(numel(fl),numel(inicios)); freq=fl; end
        area=trapz(fl,pl); if area>0,potNorm(:,k)=pl/area;else,potNorm(:,k)=pl;end
        [~,ix]=max(pl); fd(k)=refine_peak(fl,pl,ix); cl(k)=frequency_band_code(fd(k)); tm(k)=mean(t(ii))/60;
    end
    prev=nan(1,4); v=isfinite(cl); if any(v),for j=1:4,prev(j)=100*mean(cl(v)==j);end,end
    rsa=struct('tempoMin',tm,'freqCpm',freq,'potNorm',potNorm,'freqDomCpm',fd,'classe',cl,'prevalencia',prev,'janelaN',wn,'passoN',passo,'nfft',nfft);
end

function c=frequency_band_code(f)
    if f>=1&&f<2,c=1;elseif f>=2&&f<=4,c=2;elseif f>4&&f<9,c=3;elseif f>=9&&f<=12,c=4;else,c=NaN;end
end

function s=classify_frequency(f)
    c=frequency_band_code(f); names={'Bradygastria','Normogastria','Tachygastria','Intestinal band'}; if isfinite(c),s=names{c};else,s='Outside physiological bands';end
end

function fig=plot_fft_welch(fi,wi,cfg,classe,fwhm)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); subplot(2,1,1); plot(fi.fLocalCpm,fi.aLocal,'LineWidth',1.1); hold on; plot(fi.picoCpm,fi.amplitudePico,'o'); grid on; xlim([cfg.cpmMin cfg.cpmMax]); xlabel('Frequency (cpm)');ylabel('Amplitude (pT)');title(sprintf('Global FFT - peak %.3f cpm',fi.picoCpm));
    subplot(2,1,2); plot(wi.fLocalCpm,wi.pLocal,'LineWidth',1.1);hold on;plot(wi.picoCpm,wi.potenciaPico,'o'); if all(isfinite(fwhm(2:3))), yl=ylim; plot([fwhm(2) fwhm(2)],yl,'--'); plot([fwhm(3) fwhm(3)],yl,'--');end; grid on;xlim([cfg.cpmMin cfg.cpmMax]);xlabel('Frequency (cpm)');ylabel('PSD (pT^2/Hz)');title(sprintf('Global Welch PSD - peak %.3f cpm - %s',wi.picoCpm,classe));
end

function fig=plot_rsa(t,x,rsa,cfg)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); subplot(2,1,1); imagesc(rsa.tempoMin,rsa.freqCpm,rsa.potNorm); axis xy; colorbar; hold on; for y=[2 4 9], plot(xlim,[y y],'--'); end; ylim([cfg.cpmMin cfg.cpmMax]); ylabel('Frequency (cpm)'); title('Global RSA normalized per window');
    subplot(2,1,2); plot(t/60,x,'LineWidth',0.8); grid on; xlabel('Time (min)');ylabel('Amplitude (pT)');title('Global mean of valid channels');xlim([t(1) t(end)]/60);
end

function fig=plot_classification(fraction,prevalence)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); names={'Brady','Normo','Tachy','Intestinal'}; subplot(2,1,1);bar(fraction);set(gca,'XTick',1:4,'XTickLabel',names);ylabel('Power fraction (%)');title('Global power distribution');grid on; subplot(2,1,2);bar(prevalence);set(gca,'XTick',1:4,'XTickLabel',names);ylabel('RSA prevalence (%)');title('Global temporal prevalence');grid on;
end

function plot_individual_channel(ch,t,x,Fs,cfg)
    fi=calc_fft(x,Fs,cfg.cpmMin,cfg.cpmMax,cfg.nfftFFT); wi=calc_welch(x,Fs,cfg); cl=classify_frequency(wi.picoCpm); figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); subplot(3,1,1);plot(t/60,x);grid on;title(sprintf('Channel %d - time-domain signal',ch));xlabel('Time (min)');ylabel('Amplitude (pT)'); xlim([t(1) t(end)]/60); subplot(3,1,2);plot(fi.fLocalCpm,fi.aLocal);hold on;plot(fi.picoCpm,fi.amplitudePico,'o');grid on;xlim([cfg.cpmMin cfg.cpmMax]);title(sprintf('Channel %d FFT - peak %.3f cpm',ch,fi.picoCpm));xlabel('Frequency (cpm)');ylabel('Amplitude (pT)'); subplot(3,1,3);plot(wi.fLocalCpm,wi.pLocal);hold on;plot(wi.picoCpm,wi.potenciaPico,'o');grid on;xlim([cfg.cpmMin cfg.cpmMax]);title(sprintf('Channel %d Welch - peak %.3f cpm - %s',ch,wi.picoCpm,cl));xlabel('Frequency (cpm)');ylabel('PSD (pT^2/Hz)');
end

function w=hann_compat(n)
    if exist('hann','file')==2, w=hann(n,'periodic'); else, w=hanning(n); end
end

function p=percentile_compat(x,q)
    x=x(isfinite(x)); if isempty(x),p=NaN;return;end; x=sort(x); pos=1+(numel(x)-1)*q/100; lo=floor(pos);hi=ceil(pos); if lo==hi,p=x(lo);else,p=x(lo)+(pos-lo)*(x(hi)-x(lo));end
end

function s=strjoin_compatible(c,delim)
    if isempty(c),s='';return;end; s=c{1}; for i=2:numel(c),s=[s delim c{i}];end
end