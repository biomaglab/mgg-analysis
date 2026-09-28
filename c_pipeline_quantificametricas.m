function pipeline_bsgm_pre_posGari(filePre, filePos, fMinUser, fMaxUser)
% PIPELINE_BSGM_PRE_POSGARI
% Processa e compara gravações de até 10 min (pré e pós-prandial)
% ajustado para faixas de filtro e cálculo de GA-RI customizados pelo usuário,
% com suporte a ajuste opcional de amplitude e GA-RI pelo BMI.
%
% Inclui remoção de tendência linear (detrend) no sinal global para eliminar
% offsets/derivas de linha de base e artefatos de baixa frequência na PSD.

    if nargin < 2 || isempty(filePre) || isempty(filePos)
        [fPre, pPre] = uigetfile({'*.dat;*.txt;*.csv','Arquivos de dados'}, 'Selecione o arquivo PRE-PRANDIAL');
        if isequal(fPre,0), return; end
        filePre = fullfile(pPre, fPre);
        
        [fPos, pPos] = uigetfile({'*.dat;*.txt;*.csv','Arquivos de dados'}, 'Selecione o arquivo POS-PRANDIAL');
        if isequal(fPos,0), return; end
        filePos = fullfile(pPos, fPos);
    end

    %% 1. Configurações Globais e Entrada do Usuário (Filtro MGG + Correção BMI)
    cfg.Fs = 20.0;             % Frequência de amostragem padrão (Hz)
    cfg.nfft = 32768;
    cfg.fatorPicoTesla = 1.0;   % Fator unitário: o sinal de entrada já está em pT
    cfg.dpiFiguras = 300;       % Resolução visual

    % Solicita as frequências do filtro passa-banda e o valor de BMI
    if nargin < 4
        prompt = {'Qual foi a Frequência Mínima do filtro passa-banda usado no pipeline MGG? (cpm):', ...
                  'Qual foi a Frequência Máxima do filtro passa-banda usado no pipeline MGG? (cpm):', ...
                  'Informe o BMI do paciente (digite 0 ou deixe em branco para NÃO aplicar ajuste):'};
        dlgtitle = 'Parâmetros de Entrada e Ajuste Anatômico';
        dims = [1 75];
        definput = {'1.0', '6.0', '0'}; % Padrão recomendado mantido em 1.0 - 6.0 cpm
        resposta = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(resposta), return; end
        cfg.fFiltroMin = str2double(resposta{1});
        cfg.fFiltroMax = str2double(resposta{2});
        bmiInput       = str2double(resposta{3});
    else
        cfg.fFiltroMin = fMinUser;
        cfg.fFiltroMax = fMaxUser;
        
        prompt = {'Informe o BMI do paciente (digite 0 para NÃO aplicar ajuste):'};
        definput = {'0'};
        resposta = inputdlg(prompt, 'Ajuste de BMI', [1 60], definput);
        if isempty(resposta)
            bmiInput = 0;
        else
            bmiInput = str2double(resposta{1});
        end
    end

    % Trata entrada do BMI
    if isnan(bmiInput) || bmiInput <= 0
        cfg.bmi = 0;
        cfg.temBmi = false;
    else
        cfg.bmi = bmiInput;
        cfg.temBmi = true;
    end

    % Conversão das frequências limite de cpm para Hz
    cfg.fMinHz = cfg.fFiltroMin / 60.0;
    cfg.fMaxHz = cfg.fFiltroMax / 60.0;

    %% 2. Processamento dos Sinais (Pré e Pós)
    fprintf('\n---> Processando dados PRÉ-PRANDIAL (Banda: %.1f a %.1f cpm)...\n', cfg.fFiltroMin, cfg.fFiltroMax);
    resPre = processar_registro(filePre, cfg);
    
    fprintf('\n---> Processando dados PÓS-PRANDIAL (Banda: %.1f a %.1f cpm)...\n', cfg.fFiltroMin, cfg.fFiltroMax);
    resPos = processar_registro(filePos, cfg);

    %% 3. Aplicação da Correção de BMI (se informado)
    if cfg.temBmi
        f_bmi = exp(-0.05 * (cfg.bmi - 23.8));
        
        resPre.ampMedia = resPre.ampMedia / f_bmi;
        resPre.ampMax   = resPre.ampMax / f_bmi;
        resPos.ampMedia = resPos.ampMedia / f_bmi;
        resPos.ampMax   = resPos.ampMax / f_bmi;
        
        resPre.GARI = min(1.0, max(0.0, resPre.GARI / f_bmi));
        resPos.GARI = min(1.0, max(0.0, resPos.GARI / f_bmi));
        
        fprintf('[INFO] Correção por BMI (%.1f kg/m²) aplicada com sucesso.\n', cfg.bmi);
    else
        fprintf('[INFO] Processamento realizado com valores brutos (Sem ajuste de BMI).\n');
    end

    %% 4. Cálculo das Métricas Combinadas (BSGM / EGG)
    ampRatio = resPos.ampMedia / resPre.ampMedia;
    ff_AR    = resPos.ampMax / resPre.ampMedia;

    %% 5. Exibição dos Resultados no Terminal
    fprintf('\n======================================================\n');
    fprintf('           RESULTADOS DA ANÁLISE PRÉ vs PÓS           \n');
    fprintf('======================================================\n');
    fprintf('Banda do Filtro MGG : %.2f - %.2f cpm\n', cfg.fFiltroMin, cfg.fFiltroMax);
    if cfg.temBmi
        fprintf('Ajuste de BMI       : SIM (BMI = %.1f kg/m²)\n', cfg.bmi);
    else
        fprintf('Ajuste de BMI       : NÃO (Valores Brutos)\n');
    end
    fprintf('GA-RI Pré-Prandial  : %.4f (Normalidade >= 0.25)\n', resPre.GARI);
    fprintf('GA-RI Pós-Prandial  : %.4f (Normalidade >= 0.25)\n', resPos.GARI);
    fprintf('Amplitude Méd. Pré  : %.4f pT\n', resPre.ampMedia);
    fprintf('Amplitude Méd. Pós  : %.4f pT\n', resPos.ampMedia);
    fprintf('Amplitude P90 Pós   : %.4f pT\n', resPos.ampMax);
    fprintf('Amplitude Ratio     : %.4f\n', ampRatio);
    fprintf('ff-AR (Fed:Fasted)  : %.4f\n', ff_AR);
    fprintf('======================================================\n');

    %% 6. Plotagem dos Resultados
    fig = plotar_comparacao(resPre, resPos, ampRatio, ff_AR, cfg);

    %% 7. Salvar Figura e Tabela de Dados
    [folderPre, namePre, ~] = fileparts(filePre);
    outDir = fullfile(folderPre, [namePre '_comparativo_bsgm']);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    figFilePng = fullfile(outDir, 'comparativo_pre_pos_bsgm.png');
    print(fig, figFilePng, '-dpng', sprintf('-r%d', cfg.dpiFiguras));
    
    figFileMatlab = fullfile(outDir, 'comparativo_pre_pos_bsgm.fig');
    savefig(fig, figFileMatlab);

    txtFile = fullfile(outDir, 'tabela_resultados_bsgm.txt');
    fid = fopen(txtFile, 'w');
    fprintf(fid, 'METRICA\tPRE_PRANDIAL\tPOS_PRANDIAL\tVALOR_COMBINADO\tUNIDADE_STATUS\n');
    fprintf(fid, 'Banda_Filtro_cpm\t%.2f-%.2f\t%.2f-%.2f\t-\tcpm\n', cfg.fFiltroMin, cfg.fFiltroMax, cfg.fFiltroMin, cfg.fFiltroMax);
    if cfg.temBmi
        fprintf(fid, 'Ajuste_BMI\t%.1f\t%.1f\t-\tkg/m2\n', cfg.bmi, cfg.bmi);
    else
        fprintf(fid, 'Ajuste_BMI\tNao_Aplicado\tNao_Aplicado\t-\t-\n');
    end
    fprintf(fid, 'GA-RI\t%.6f\t%.6f\t-\tPre: %s | Pos: %s\n', resPre.GARI, resPos.GARI, status_gari(resPre.GARI), status_gari(resPos.GARI));
    fprintf(fid, 'Amplitude_Media\t%.6f\t%.6f\t-\tpT\n', resPre.ampMedia, resPos.ampMedia);
    fprintf(fid, 'Amplitude_P90\t%.6f\t%.6f\t-\tpT\n', resPos.ampMax, resPre.ampMedia);
    fprintf(fid, 'Amplitude_Ratio\t-\t-\t%.6f\tratio\n', ampRatio);
    fprintf(fid, 'ff-AR\t-\t-\t%.6f\tratio\n', ff_AR);
    fclose(fid);
    fprintf('[OK] Resultados e figuras salvos em: %s\n\n', outDir);
end

%% ================= FUNÇÕES AUXILIARES =================

function res = processar_registro(filePath, cfg)
    M = read_numeric_compatible(filePath);
    M = M(~all(isnan(M),2),:);
    M = M(:,~all(isnan(M),1));
    
    t = M(:,1);
    X = M(:,2:end) * cfg.fatorPicoTesla; % Mantém os valores originais em pT
    [nAmostras, nCanais] = size(X);
    
    dt = diff(t); dt = dt(isfinite(dt) & dt > 0);
    if isempty(dt), Fs = cfg.Fs; else, Fs = 1/median(dt); end
    
    Xtrat = zeros(nAmostras, nCanais);
    validos = true(1, nCanais);
    for ch = 1:nCanais
        x = X(:,ch);
        if mean(isfinite(x)) < 0.9 || max(x)-min(x) <= eps
            validos(ch) = false;
            continue;
        end
        Xtrat(:,ch) = x - mean(x);
    end
    
    canaisValidos = find(validos);
    sinalGlobal = mean(Xtrat(:, canaisValidos), 2);
    
    % =========================================================================
    % APLICAÇÃO DA OPÇÃO 1: Remoção da tendência linear (Detrending)
    % Remove a deriva/offset ascendente da linha de base
    % =========================================================================
    sinalFiltrado = detrend(sinalGlobal);

    tamJanelaSeg = 60;
    janelaWelch  = round(tamJanelaSeg * Fs);
    nOverlap     = round((tamJanelaSeg / 2) * Fs);

    [pxx, fHz] = pwelch(sinalFiltrado, hann_compat(janelaWelch), nOverlap, cfg.nfft, Fs);
    fCpm = fHz * 60;
    
    idxNarrow = (fCpm >= cfg.fFiltroMin & fCpm <= cfg.fFiltroMax);
    
    fTotalMin = min(1.0, cfg.fFiltroMin);
    fTotalMax = max(6.0, cfg.fFiltroMax);
    idxTotal  = (fCpm >= fTotalMin & fCpm <= fTotalMax);
    
    gari = sum(pxx(idxNarrow)) / max(sum(pxx(idxTotal)), eps);
    
    janelaAmp = max(2, round(20 * Fs));
    nJan = floor(nAmostras / janelaAmp);
    ampJan = zeros(nJan, 1);
    for k = 1:nJan
        seg = sinalFiltrado((k-1)*janelaAmp+1 : k*janelaAmp);
        ampJan(k) = max(seg) - min(seg);
    end
    
    res.t = t;
    res.sinal = sinalFiltrado;
    res.fCpm = fCpm;
    res.pxx = pxx;
    res.GARI = gari;
    res.ampMedia = mean(ampJan);
    res.ampMax   = prctile(ampJan, 90);
end

function fig = plotar_comparacao(resPre, resPos, ampRatio, ff_AR, cfg)
    fig = figure('Color', 'w', 'Name', 'Comparativo BSGM PRE vs POS', 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
    
    minGlobal = min([resPre.sinal; resPos.sinal]);
    maxGlobal = max([resPre.sinal; resPos.sinal]);
    margem = 0.05 * (maxGlobal - minGlobal);
    if margem == 0, margem = 1; end
    limitesY = [minGlobal - margem, maxGlobal + margem];

    subplot(2,2,1);
    plot(resPre.t / 60, resPre.sinal, 'b', 'LineWidth', 1);
    title('Sinal Global Pré-Prandial');
    xlabel('Tempo (min)'); ylabel('Amplitude (pT)'); grid on;
    ylim(limitesY);
    
    subplot(2,2,2);
    plot(resPos.t / 60, resPos.sinal, 'r', 'LineWidth', 1);
    title('Sinal Global Pós-Prandial');
    xlabel('Tempo (min)'); ylabel('Amplitude (pT)'); grid on;
    ylim(limitesY);
    
    % Subplot 3: Espectro de Potência (PSD)
    subplot(2,2,3);
    hold on; grid on;
    
    plot(resPre.fCpm, resPre.pxx, 'b', 'LineWidth', 1.5, 'DisplayName', 'Pré-Prandial');
    plot(resPos.fCpm, resPos.pxx, 'r', 'LineWidth', 1.5, 'DisplayName', 'Pós-Prandial');
    
    % Restrição do pico para a janela de Normogastria (2.0 - 4.0 cpm)
    idxNormo = (resPre.fCpm >= 2.0 & resPre.fCpm <= 4.0);
    fNormo = resPre.fCpm(idxNormo);
    
    if ~isempty(fNormo)
        [maxPre, idxPre] = max(resPre.pxx(idxNormo));
        fPicoPre = fNormo(idxPre);
        
        [maxPos, idxPos] = max(resPos.pxx(idxNormo));
        fPicoPos = fNormo(idxPos);
        
        plot(fPicoPre, maxPre, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(fPicoPos, maxPos, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
        
        text(fPicoPre, maxPre * 1.04, sprintf('  Pré Pico: %.2f cpm', fPicoPre), 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 8.5);
        text(fPicoPos, maxPos * 1.04, sprintf('  Pós Pico: %.2f cpm', fPicoPos), 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 8.5);
    end
        
    xlim([min(1.0, cfg.fFiltroMin - 0.2), max(6.0, cfg.fFiltroMax + 0.5)]);
    
    idxBandaFiltro = (resPre.fCpm >= cfg.fFiltroMin & resPre.fCpm <= cfg.fFiltroMax);
    maxLocal = max([resPre.pxx(idxBandaFiltro); resPos.pxx(idxBandaFiltro)]);
    if maxLocal > 0, ylim([0, maxLocal * 1.25]); end
    
    title(sprintf('PSD (Filtro: %.1f - %.1f cpm | Pico em Normogastria)', cfg.fFiltroMin, cfg.fFiltroMax));
    xlabel('Frequência (cpm)'); ylabel('PSD (pT^2/Hz)');
    legend('Location', 'northeast');
    
    subplot(2,2,4);
    axis off;
    
    txtPre = status_gari(resPre.GARI);
    txtPos = status_gari(resPos.GARI);
    
    if cfg.temBmi
        txtBmi = sprintf('Ajuste BMI: SIM (%.1f kg/m²)', cfg.bmi);
    else
        txtBmi = 'Ajuste BMI: NÃO (Bruto)';
    end
    
    strText = {
        'Métricas de BSGM e EGG:';
        sprintf('Banda Filtro MGG: %.2f - %.2f cpm', cfg.fFiltroMin, cfg.fFiltroMax);
        txtBmi;
        '';
        sprintf('• GA-RI Pré: %.3f %s', resPre.GARI, txtPre);
        sprintf('• GA-RI Pós: %.3f %s', resPos.GARI, txtPos);
        '';
        sprintf('• Amplitude Méd. Pré: %.2f pT', resPre.ampMedia);
        sprintf('• Amplitude Méd. Pós: %.2f pT', resPos.ampMedia);
        sprintf('• Amplitude P90 Pós: %.2f pT', resPos.ampMax);
        '';
        sprintf('• Amplitude Ratio: %.2f', ampRatio);
        sprintf('• ff-AR (Fed:Fasted P90): %.2f', ff_AR)
    };
    
    text(0.05, 0.5, strText, 'Interpreter', 'none', 'FontSize', 10);
end

function txt = status_gari(val)
    if val >= 0.25
        txt = '(Normal >= 0.25)';
    else
        txt = '(Anormal < 0.25)';
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