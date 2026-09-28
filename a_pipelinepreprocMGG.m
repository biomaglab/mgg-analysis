function pipeline_EGG_tutorial_sem_ICA_inspecao_iterativa(inputFile, Fs)
% PIPELINE_EGG_TUTORIAL_SEM_ICA_INSPECAO_ITERATIVA
% Pipeline multicanal para arquivos .dat de EGG/MGG com amplitude em picoTesla (pT).

    if nargin < 1 || isempty(inputFile)
        [file, folder] = uigetfile( ...
            {'*.dat;*.txt;*.csv;*.tsv', 'Dados numericos'; '*.*', 'Todos'}, ...
            'Selecione o arquivo .dat');
        if isequal(file, 0)
            fprintf('Nenhum arquivo selecionado.\n');
            return;
        end
        inputFile = fullfile(folder, file);
    end

    %% Configuracao
    cfg.firstColumnIsTime = true;
    cfg.defaultFs = 20;
    cfg.nSigma = 6.0;                  % limiar conservador: 6 MAD
    cfg.hampelWindowSeconds = 0.55;   % janela total alvo do Hampel (~11 amostras em 20 Hz)
    cfg.hampelHalfWindow = [];         % vazio = calculado automaticamente a partir de hampelWindowSeconds
    cfg.defaultLowHz = 0.015;
    cfg.defaultHighHz = 0.2;
    cfg.butterOrder = 5;
    cfg.movementFrequencyCpm = 3.0;
    cfg.movementWindowCycles = 1.0;
    cfg.fatorPicoTesla = 400.0;        % Fator de conversao para picoTesla (pT) [*400]

    %% Carregamento
    M = read_numeric_compatible(inputFile);
    
    M = M(~all(isnan(M), 2), :);
    M = M(:, ~all(isnan(M), 1));

    if isempty(M) || size(M,1) < 2
        error('O arquivo nao contem uma matriz numerica valida.');
    end

    if cfg.firstColumnIsTime && size(M,2) >= 2
        t = M(:,1);
        X = M(:,2:end) * cfg.fatorPicoTesla; % Conversao para pT aplicada (*400)
        dt = diff(t);
        dt = dt(isfinite(dt) & dt > 0);

        if nargin < 2 || isempty(Fs)
            if isempty(dt)
                Fs = cfg.defaultFs;
                warning('Nao foi possivel inferir Fs pelo tempo. Usando %.6g Hz.', Fs);
            else
                Fs = 1 / median(dt);
            end
        end
    else
        X = M * cfg.fatorPicoTesla; % Conversao para pT aplicada (*400)
        if nargin < 2 || isempty(Fs)
            Fs = cfg.defaultFs;
        end
        t = (0:size(X,1)-1)' / Fs;
    end

    if ~isfinite(Fs) || Fs <= 0
        error('Fs deve ser positivo.');
    end

    fprintf('\nArquivo: %s\n', inputFile);
    fprintf('Amostras: %d | canais: %d | Fs: %.8g Hz | duracao: %.2f min\n', ...
        size(X,1), size(X,2), Fs, (t(end)-t(1))/60);

    %% Corrige apenas NaN/Inf
    for ch = 1:size(X,2)
        good = isfinite(X(:,ch));
        if nnz(good) < 2
            error('Canal %d contem menos de duas amostras validas.', ch);
        end
        if ~all(good)
            X(~good,ch) = interp1(find(good), X(good,ch), ...
                find(~good), 'linear', 'extrap');
        end
    end
    raw = X;

    %% Escolhas iniciais do usuario
    cfg.outlierMethod = solicitar_metodo_outlier();
    [cfg.lowHz, cfg.highHz] = solicitar_banda_usuario( ...
        Fs, cfg.defaultLowHz, cfg.defaultHighHz);

    %% 1. Centralizacao robusta pela mediana
    centered = bsxfun(@minus, raw, nanmedian_compatible(raw,1));
    movementHz = cfg.movementFrequencyCpm / 60;

    %% Processamento iterativo com inspecao ANTES de salvar
    aceito = false;
    tentativa = 0;
    figsInspecao = [];

    while ~aceito
        tentativa = tentativa + 1;

        if ~isempty(figsInspecao)
            for kk = 1:numel(figsInspecao)
                if ishandle(figsInspecao(kk))
                    close(figsInspecao(kk));
                end
            end
        end

        fprintf('\n============================================================\n');
        fprintf('PROCESSAMENTO / INSPECAO - TENTATIVA %d\n', tentativa);
        fprintf('============================================================\n');
        fprintf('Outliers: %s | nSigma = %.3g', upper(cfg.outlierMethod), cfg.nSigma);
        if strcmpi(cfg.outlierMethod,'hampel')
            fprintf(' | janela alvo = %.3g s', cfg.hampelWindowSeconds);
        end
        fprintf('\n');
        fprintf('Butterworth: %.6g-%.6g Hz (%.4g-%.4g cpm), ordem %d, zero-phase\n', ...
            cfg.lowHz, cfg.highHz, cfg.lowHz*60, cfg.highHz*60, cfg.butterOrder);

        %% 2. Hampel ou MAD
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
                fprintf('Hampel: janela = %d amostras (%.3f s), meia janela k=%d, limiar=%.3g MAD\n', ...
                    nWin, nWin/Fs, k, cfg.nSigma);
                [despiked, hampelCount] = hampel_multichannel(centered, k, cfg.nSigma);
                fprintf('Hampel corrigiu %d de %d pontos (%.3f%%).\n', ...
                    hampelCount, numel(centered), 100*hampelCount/numel(centered));
            case 'mad'
                despiked = mad_filter_multichannel(centered, cfg.nSigma);
            otherwise
                error('Metodo de outlier invalido.');
        end

        %% 3. Butterworth definido pelo usuario
        lowHz = cfg.lowHz;
        highHz = cfg.highHz;
        if ~(0 < lowHz && lowHz < highHz && highHz < Fs/2)
            error('Banda %.6g-%.6g Hz invalida para Fs=%.6g Hz.', lowHz, highHz, Fs);
        end

        [z,p,kFilt] = butter(cfg.butterOrder, [lowHz highHz]/(Fs/2), 'bandpass');
        sos = zp2sos(z,p,kFilt);
        bandpassed = zeros(size(despiked));
        for ch = 1:size(despiked,2)
            bandpassed(:,ch) = sosfiltfilt_compatible(sos, despiked(:,ch));
        end

        %% 4. LMMSE de movimento
        [movementClean, movementNoise] = movement_filter_lmmse( ...
            bandpassed, Fs, movementHz, cfg.movementWindowCycles);

        %% Inspecao VISIVEL
        figsInspecao = create_figures_inspection(t, raw, centered, despiked, ...
            bandpassed, movementClean, cfg, tentativa);

        fprintf('\nAs figuras estao abertas para inspecao.\n');
        fprintf('Nenhum arquivo de resultado foi salvo ainda.\n');

        escolha = solicitar_decisao_inspecao(cfg.outlierMethod);

        switch escolha
            case 1
                aceito = true;
            case 2
                if strcmpi(cfg.outlierMethod,'hampel')
                    cfg.outlierMethod = 'mad';
                else
                    cfg.outlierMethod = 'hampel';
                end
                fprintf('\nMetodo alterado para %s. Reprocessando...\n', upper(cfg.outlierMethod));
            case 3
                [cfg.lowHz, cfg.highHz] = solicitar_banda_usuario(Fs, cfg.lowHz, cfg.highHz);
                fprintf('\nBanda alterada. Reprocessando...\n');
            case 4
                fprintf('\nProcessamento cancelado pelo usuario.\n');
                if ~isempty(figsInspecao)
                    for kk = 1:numel(figsInspecao)
                        if ishandle(figsInspecao(kk)), close(figsInspecao(kk)); end
                    end
                end
                return;
        end
    end

    %% Saida - NOMENCLATURA E SALVAMENTO DE ARQUIVOS
    [folder, base, ~] = fileparts(inputFile);

    % Limpa o nome removendo os trechos indesejados
    baseLimpa = strrep(base, '_original_dec20Hz', '');
    baseLimpa = strrep(baseLimpa, '_pipeline_sem_ICA', '');

    % Verifica se contem 'pos' usando strfind (compativel com todas as versoes do MATLAB)
    if ~isempty(strfind(lower(baseLimpa), 'pos'))
        sufixo = '(pos)';
    else
        sufixo = '(pre)';
    end

    outputDir = fullfile(folder, [baseLimpa sufixo]);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % SALVAMENTO DOS ARQUIVOS NA PASTA
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
    save_inspection_figures(figsInspecao, outputDir);

    fprintf('\nPipeline concluido com sucesso.\n');
    fprintf('Resultados salvos em: %s\n', outputDir);
end

%% ================= FUNÇÕES AUXILIARES =================

function method = solicitar_metodo_outlier()
    fprintf('\nRejeicao de outliers:\n');
    fprintf('  1 - Hampel local\n');
    fprintf('  2 - MAD global\n');
    s = input('Escolha [Enter = 1]: ', 's');
    if isempty(strtrim(s)) || strcmp(strtrim(s),'1')
        method = 'hampel';
    elseif strcmp(strtrim(s),'2')
        method = 'mad';
    else
        error('Opcao invalida. Use 1 ou 2.');
    end
end

function [lowHz, highHz] = solicitar_banda_usuario(Fs, defaultLowHz, defaultHighHz)
    fprintf('\nDefina o Butterworth passa-faixa em Hz.\n');
    s1 = input(sprintf('Frequencia inferior [%.6g Hz]: ', defaultLowHz), 's');
    s2 = input(sprintf('Frequencia superior [%.6g Hz]: ', defaultHighHz), 's');

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
        error('As frequencias devem ser numericas.');
    end
    if ~(lowHz > 0 && lowHz < highHz && highHz < Fs/2)
        error('Banda invalida.');
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
        error('Frequencia e comprimento da janela devem ser positivos.');
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
    fprintf(fid,'PIPELINE EGG/MGG - SEM ICA\n');
    fprintf(fid,'Arquivo: %s\n',inputFile);
    fprintf(fid,'Fs: %.10g Hz\n',Fs);
    fprintf(fid,'Canais: %d\n',size(raw,2));
    fprintf(fid,'Metodo de outliers: %s\n',upper(cfg.outlierMethod));
    fprintf(fid,'Butterworth: %.6g-%.6g Hz (%.4g-%.4g cpm)\n', cfg.lowHz,cfg.highHz,cfg.lowHz*60,cfg.highHz*60);
    fprintf(fid,'Sinal final: 05_movement_filtered.dat\n');
end

function escolha = solicitar_decisao_inspecao(method)
    fprintf('\nDECISAO APOS INSPECAO DOS GRAFICOS:\n');
    fprintf('  1 - Aceitar e SALVAR\n');
    fprintf('  2 - Trocar metodo de outlier e reavaliar\n');
    fprintf('  3 - Alterar banda Butterworth e reavaliar\n');
    fprintf('  4 - Cancelar sem salvar\n');

    while true
        s = input('Escolha [1-4]: ','s');
        escolha = str2double(strtrim(s));
        if isfinite(escolha) && any(escolha == 1:4)
            return;
        end
        fprintf('Opcao invalida. Digite 1, 2, 3 ou 4.\n');
    end
end

function figs = create_figures_inspection(t, raw, centered, despiked, bandpassed, movementClean, cfg, tentativa)
    stages = { ...
        'Sinal original', raw; ...
        'Apos centralizacao pela mediana', centered; ...
        sprintf('Apos %s', upper(cfg.outlierMethod)), despiked; ...
        sprintf('Butterworth %.4g-%.4g Hz',cfg.lowHz,cfg.highHz), bandpassed; ...
        'Apos LMMSE', movementClean};

    figs = zeros(1,size(stages,1));
    for s = 1:size(stages,1)
        X = stages{s,2};
        nCh = size(X,2);
        figs(s) = figure_fullscreen_compat('Visible','on','Color','w', ...
            'Name',sprintf('Tentativa %d - %s',tentativa,stages{s,1}), ...
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
                xlabel(ax,'Tempo (min)');
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
    nomes = { ...
        'figure_01_raw.png', ...
        'figure_02_median_centered.png', ...
        'figure_03_outlier_filtered.png', ...
        'figure_04_bandpass.png', ...
        'figure_05_movement_filtered.png'};

    for k = 1:min(numel(figs), numel(nomes))
        if isgraphics(figs(k))
            print(figs(k), fullfile(outputDir, nomes{k}), '-dpng', '-r160');
        end
    end
end

function y = sosfiltfilt_compatible(sos, x)
    x = x(:);
    y = x;
    for sec = 1:size(sos,1)
        b = sos(sec,1:3);
        a = sos(sec,4:6);
        if abs(a(1)) < eps, error('Secao SOS invalida.'); end
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
            error('O arquivo nao contem matriz numerica.');
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