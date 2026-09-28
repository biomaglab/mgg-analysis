function analiseglobal_pospipeline_inspecao_canais_fft_welch(inputFile)
% ANALISEGLOBAL_POSPIPELINE_INSPECAO_CANAIS_FFT_WELCH
% Processamento e análise global de sinais bioelétricos/biomagnéticos gástricos.
%
% Recursos incluídos:
% - Leitura direta na escala original em picoTesla (pT).
% - Seleção de corte de tempo inicial via cursor/cruzeta interativo (ginput).
% - Controle de qualidade preliminar por canal (dados, platô, sinal constante, RMS).
% - Análise espectral global por FFT, Welch e RSA (Running Spectral Analysis).
% - Cálculo do Instability Coefficient (IC) baseado no modelo GastroPy.
% - Exportação automática dos dados em arquivos .DAT, .TXT, .FIG e .PNG.

    if nargin < 1 || isempty(inputFile)
        [file, folder] = uigetfile({'*.dat;*.txt;*.csv;*.tsv','Arquivos numéricos'; '*.*','Todos'}, ...
            'Selecione o arquivo pós-pipeline');
        if isequal(file,0)
            fprintf('Nenhum arquivo selecionado.\n');
            return;
        end
        inputFile = fullfile(folder,file);
    end

    %% Configuração do Pipeline
    cfg.fsPadrao = 20.0;
    cfg.freqMinHz = 0.015;
    cfg.freqMaxHz = 0.200;
    cfg.tInicioSeg = 50.0; % Tempo inicial padrão descartado
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
    cfg.fatorPicoTesla = 1.0; % Fator unitário: o sinal já está em pT

    BRADI = [1 2]; NORMO = [2 4]; TAQUI = [4 9]; INTESTINAL = [9 12];
    nomesFaixas = {'Bradigastria','Normogastria','Taquigastria','Banda intestinal'};

    %% Leitura dos Dados
    M = read_numeric_compatible(inputFile);
    M = M(~all(isnan(M),2),:);
    M = M(:,~all(isnan(M),1));
    if size(M,2) < 2 || size(M,1) < 2
        error('Arquivo inválido: esperado vetor de tempo na coluna 1 e pelo menos um canal.');
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
    if ~isfinite(Fs) || Fs <= 0, error('Frequência de amostragem inválida.'); end

    fprintf('\nArquivo: %s\n', inputFile);
    fprintf('Fs detectada: %.8g Hz\n', Fs);
    fprintf('Duração total original: %.3f s (%.3f min)\n', tempoOrig(end)-tempoOrig(1), (tempoOrig(end)-tempoOrig(1))/60);

    %% 1. Inspeção Prévia e Corte Temporal por Cursor (Cruzeta)
    [tClicado, figPrevia] = plot_sinal_bruto_previa_com_cruzeta(tempoOrig, XorigBruto, cfg.tInicioSeg);
    
    if isfinite(tClicado)
        cfg.tInicioSeg = tClicado;
    end

    fprintf('\n-----------------------------------------------------------------\n');
    fprintf('Tempo selecionado visualmente: %.2f s\n', cfg.tInicioSeg);
    fprintf('-----------------------------------------------------------------\n');

    r = input(sprintf('Pressione Enter para confirmar %.2f s ou digite um novo valor (0 = sem corte): ', cfg.tInicioSeg), 's');
    if ~isempty(strtrim(r))
        valInput = str2double(r);
        if isfinite(valInput) && valInput >= 0
            cfg.tInicioSeg = valInput;
        end
    end

    if isgraphics(figPrevia) && isvalid(figPrevia)
        close(figPrevia);
    end

    % Aplicação do corte de tempo inicial
    tempo = tempoOrig;
    Xorig = XorigBruto;
    if cfg.tInicioSeg > 0
        idxInicio = find(tempoOrig >= tempoOrig(1) + cfg.tInicioSeg, 1, 'first');
        if ~isempty(idxInicio) && idxInicio < numel(tempoOrig)
            tempo = tempoOrig(idxInicio:end);
            Xorig = XorigBruto(idxInicio:end, :);
            fprintf('--> Sinal cortado! Analisando de t = %.2f s (%.2f min) até t = %.2f s.\n', tempo(1), tempo(1)/60, tempo(end));
        else
            warning('Tempo de corte informado é maior que a duração do sinal. O corte foi ignorado.');
        end
    else
        fprintf('--> Nenhum corte aplicado. Analisando todo o sinal a partir de t = %.2f s.\n', tempo(1));
    end

    [nAmostras, nCanais] = size(Xorig);
    fprintf('Amostras úteis: %d | Canais: %d\n', nAmostras, nCanais);
    fprintf('Duração útil: %.3f s (%.3f min)\n', tempo(end)-tempo(1), (tempo(end)-tempo(1))/60);
    fprintf('Escala de amplitude: Mantida na escala de entrada (picoTesla - pT)\n');

    %% Definir Banda Espectral da Análise
    fprintf('\nA rotina NÃO reaplica filtro passa-faixa.\n');
    fprintf('Informe apenas a faixa espectral que deseja analisar.\n');
    r = input(sprintf('Frequência inferior [%.6g Hz = %.3f cpm]: ', cfg.freqMinHz, cfg.freqMinHz*60), 's');
    if ~isempty(strtrim(r)), cfg.freqMinHz = str2double(r); end
    r = input(sprintf('Frequência superior [%.6g Hz = %.3f cpm]: ', cfg.freqMaxHz, cfg.freqMaxHz*60), 's');
    if ~isempty(strtrim(r)), cfg.freqMaxHz = str2double(r); end
    if ~(isfinite(cfg.freqMinHz)&&isfinite(cfg.freqMaxHz)&&0<cfg.freqMinHz&&cfg.freqMinHz<cfg.freqMaxHz&&cfg.freqMaxHz<Fs/2)
        error('Faixa espectral inválida.');
    end
    cfg.cpmMin = cfg.freqMinHz*60;
    cfg.cpmMax = cfg.freqMaxHz*60;
    fprintf('Faixa analisada: %.4f-%.4f cpm (%.6g-%.6g Hz).\n', cfg.cpmMin, cfg.cpmMax, cfg.freqMinHz, cfg.freqMaxHz);

    %% Controle Técnico por Canal
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
        maiorSeq=maior_sequencia_true(repetido);
        maiorPlatoSeg(ch)=(maiorSeq+1)/Fs;
        falhaPlato(ch)=maiorPlatoSeg(ch)>=cfg.duracaoPlatoMaxSeg;
        Xtrat(:,ch)=x;
        rmsCanal(ch)=sqrt(mean(x.^2));
    end

    candidatos=find(~falhaDados & ~falhaConstante & ~falhaPlato & isfinite(rmsCanal));
    if isempty(candidatos), error('Nenhum canal tecnicamente válido foi encontrado.'); end
    medRMS=median(rmsCanal(candidatos));
    limiarRMS=cfg.fatorRMSQuaseNulo*medRMS;
    falhaRMS(candidatos)=rmsCanal(candidatos)<limiarRMS;
    alertaRMS=falhaRMS;
    falhaAuto=falhaDados|falhaConstante|falhaPlato;

    fprintf('\nControle técnico preliminar:\n');
    for ch=1:nCanais
        avisos={};
        if falhaDados(ch), avisos{end+1}='dados insuficientes'; end
        if falhaConstante(ch), avisos{end+1}='sinal constante'; end
        if falhaPlato(ch), avisos{end+1}='platô prolongado'; end
        if alertaRMS(ch), avisos{end+1}='RMS quase nulo - ALERTA'; end
        if isempty(avisos), txt='sem alerta técnico'; else txt=strjoin_compatible(avisos,'; '); end
        fprintf('  Canal %d: %s\n',ch,txt);
    end

    %% Inspeção Visual dos Canais
    figInspecao=plot_inspecao(tempo,Xorig,falhaDados,falhaConstante,falhaPlato,alertaRMS);
    fprintf('\nInspecione a figura de qualidade dos canais. Feche a janela para continuar.\n');
    waitfor(figInspecao);

    fprintf('\nCanais com falha objetiva (dados/constante/platô) serão excluídos automaticamente.\n');
    fprintf('RMS baixo é apenas alerta e NÃO exclui automaticamente.\n');
    resposta=input('Canais adicionais a excluir [ex: 2 5 7 | Enter = nenhum]: ','s');
    exclusaoManual=parse_channel_list(resposta,nCanais);

    canaisValidosIdx=find(~falhaAuto & ~ismember(1:nCanais,exclusaoManual));
    canaisExcluidosIdx=setdiff(1:nCanais,canaisValidosIdx);
    if isempty(canaisValidosIdx), error('Nenhum canal válido permaneceu para a análise global.'); end
    fprintf('\nCanais usados na média global: %s\n',num2str(canaisValidosIdx));
    if isempty(canaisExcluidosIdx), fprintf('Canais excluídos: nenhum\n'); else fprintf('Canais excluídos: %s\n',num2str(canaisExcluidosIdx)); end

    %% Sinal Global Médio e Mediana
    sinalGlobal=mean(Xtrat(:,canaisValidosIdx),2);
    sinalGlobal=sinalGlobal-mean(sinalGlobal);
    sinalGlobalMediana=median(Xtrat(:,canaisValidosIdx),2);
    sinalGlobalMediana=sinalGlobalMediana-mean(sinalGlobalMediana);

    %% FFT + Welch Global
    fftInfo=calc_fft(sinalGlobal,Fs,cfg.cpmMin,cfg.cpmMax,cfg.nfftFFT);
    welchInfo=calc_welch(sinalGlobal,Fs,cfg);
    classeGlobal=classificar_freq(welchInfo.picoCpm);

    potenciaTotal=trapz(welchInfo.fLocalCpm,welchInfo.pLocal);
    potFaixas=[integrar_banda(welchInfo.fCpm,welchInfo.pxx,BRADI), ...
               integrar_banda(welchInfo.fCpm,welchInfo.pxx,NORMO), ...
               integrar_banda(welchInfo.fCpm,welchInfo.pxx,TAQUI), ...
               integrar_banda(welchInfo.fCpm,welchInfo.pxx,INTESTINAL)];
    if sum(potFaixas)>0, fracaoPotencia=100*potFaixas/sum(potFaixas); else, fracaoPotencia=nan(1,4); end
    fwhm=calc_fwhm(welchInfo.fLocalCpm,welchInfo.pLocal,welchInfo.idxPicoLocal);
    entropia=entropia_espectral(welchInfo.pLocal);

    %% Cálculo do Instability Coefficient (IC)
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

    %% Controle Individual por Canal
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
            controle(ch).classe_fft=classificar_freq(fi.picoCpm);
            controle(ch).classe_welch=classificar_freq(wi.picoCpm);
            controle(ch).delta=abs(fi.picoCpm-wi.picoCpm);
            if controle(ch).delta<=limiarConcordancia
                controle(ch).concordancia='CONCORDAM';
            elseif strcmp(controle(ch).classe_fft,controle(ch).classe_welch)
                controle(ch).concordancia='MESMA_FAIXA_PICOS_DIVERGENTES';
            else
                controle(ch).concordancia='DIVERGEM_ENTRE_FAIXAS';
            end
            controle(ch).entropia=entropia_espectral(wi.pLocal);
            controle(ch).ic=icChInfo.IC;
            if alertaRMS(ch), controle(ch).status='VALIDO_COM_ALERTA_RMS'; else, controle(ch).status='VALIDO'; end
        else
            controle(ch).fft_cpm=NaN; controle(ch).welch_cpm=NaN;
            controle(ch).classe_fft='NAO_ANALISADO'; controle(ch).classe_welch='NAO_ANALISADO';
            controle(ch).delta=NaN; controle(ch).concordancia='NAO_ANALISADO'; controle(ch).entropia=NaN;
            controle(ch).ic=NaN;
            if ismember(ch,exclusaoManual), controle(ch).status='EXCLUSAO_MANUAL_FALHA_CONFIRMADA';
            elseif falhaDados(ch), controle(ch).status='FALHA_DADOS';
            elseif falhaConstante(ch), controle(ch).status='SINAL_CONSTANTE';
            elseif falhaPlato(ch), controle(ch).status='PLATO_PROLONGADO';
            else, controle(ch).status='EXCLUIDO'; end
        end
    end

    %% Diretório de Saída
    [folder, base, ~] = fileparts(inputFile);
    baseLimpa = strrep(base, '_original_dec20Hz', '');
    baseLimpa = strrep(baseLimpa, '_pipeline_sem_ICA', '');
    baseLimpa = strrep(baseLimpa, '_analise_global', '');

    if ~isempty(strfind(lower(baseLimpa), 'pos'))
        sufixo = '(pos)';
    else
        sufixo = '(pre)';
    end

    outDir = fullfile(folder, [baseLimpa sufixo]);
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    %% Geração de Figuras
    figs(1) = plot_empilhados(tempo, Xtrat, canaisValidosIdx);
    figs(2) = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
    plot(tempo/60, sinalGlobal, 'LineWidth', 0.9); grid on; 
    xlabel('Tempo (min)'); ylabel('Amplitude (pT)'); 
    title(sprintf('Traçado global médio - %d canais válidos (Iniciando em %.1fs)', numel(canaisValidosIdx), tempo(1)));

    figs(3) = plot_fft_welch(fftInfo, welchInfo, cfg, classeGlobal, fwhm);
    figs(4) = plot_rsa(tempo, sinalGlobal, rsa, cfg);
    figs(5) = plot_classificacao(fracaoPotencia, rsa.prevalencia);

    nomesBase = {'01_sinais_bandpass', ...
                 '02_tracado_global', ...
                 '03_fft_welch_global', ...
                 '04_rsa_global_e_tracado', ...
                 '05_classificacao_fisiologica'};

    for k = 1:numel(figs)
        if isgraphics(figs(k)) && isvalid(figs(k))
            print(figs(k), fullfile(outDir, [nomesBase{k} '.png']), '-dpng', sprintf('-r%d', cfg.dpiFiguras));
            savefig(figs(k), fullfile(outDir, [nomesBase{k} '.fig']));
        end
    end

    %% Exportação de Dados numéricos em arquivos .DAT
    % A) Traçado Temporal Global em .DAT
    matrizSinaisDAT = [tempo, tempo/60, sinalGlobal, sinalGlobalMediana];
    fileSinaisDAT = fullfile(outDir, 'sinal_global_tempo.dat');
    fidDAT = fopen(fileSinaisDAT, 'w');
    fprintf(fidDAT, 'Tempo_s\tTempo_min\tSinal_Media_pT\tSinal_Mediana_pT\n');
    fclose(fidDAT);
    dlmwrite(fileSinaisDAT, matrizSinaisDAT, '-append', 'delimiter', '\t', 'precision', '%.8f');

    % B) Espectro de Frequência Welch e FFT em .DAT
    matrizEspectroDAT = [welchInfo.fLocalCpm/60, welchInfo.fLocalCpm, welchInfo.pLocal, fftInfo.aLocal];
    fileEspectroDAT = fullfile(outDir, 'espectro_welch_fft.dat');
    fidDAT = fopen(fileEspectroDAT, 'w');
    fprintf(fidDAT, 'Freq_Hz\tFreq_cpm\tPSD_Welch_pT2_Hz\tAmp_FFT_pT\n');
    fclose(fidDAT);
    dlmwrite(fileEspectroDAT, matrizEspectroDAT, '-append', 'delimiter', '\t', 'precision', '%.8f');

    % C) Tabela de Controle dos Canais em .DAT
    fileControleDAT = fullfile(outDir, 'controle_canais.dat');
    fidDAT = fopen(fileControleDAT, 'w');
    fprintf(fidDAT, 'Canal\tValido_pct\tRMS_pT\tFFT_cpm\tClasse_FFT\tWelch_cpm\tClasse_Welch\tDelta\tConcordancia\tEntropia\tIC\tStatus\n');
    for ch=1:nCanais
        fprintf(fidDAT, '%d\t%.6f\t%.10g\t%.8f\t%s\t%.8f\t%s\t%.8f\t%s\t%.8f\t%.8f\t%s\n', ...
            controle(ch).canal, controle(ch).valido_pct, controle(ch).rms, controle(ch).fft_cpm, controle(ch).classe_fft, ...
            controle(ch).welch_cpm, controle(ch).classe_welch, controle(ch).delta, controle(ch).concordancia, controle(ch).entropia, controle(ch).ic, controle(ch).status);
    end
    fclose(fidDAT);

    % D) Sinais Filtrados dos 10 Canais em .DAT
    matrizCanaisDAT = [tempo, Xtrat];
    fileCanaisDAT = fullfile(outDir, 'sinais_canais_filtrados.dat');
    fidDAT = fopen(fileCanaisDAT, 'w');
    headersCanais = sprintf('Tempo_s\t');
    for ch = 1:nCanais
        headersCanais = [headersCanais, sprintf('Canal_%d_pT\t', ch)]; %#ok<AGROW>
    end
    headersCanais = [strtrim(headersCanais), '\n'];
    fprintf(fidDAT, headersCanais);
    fclose(fidDAT);
    dlmwrite(fileCanaisDAT, matrizCanaisDAT, '-append', 'delimiter', '\t', 'precision', '%.8f');

    %% Relatório TXT
    fid=fopen(fullfile(outDir,'relatorio_analise_global.txt'),'w');
    fprintf(fid,'ANALISE GLOBAL GASTRICA POS-PIPELINE\n===================================\n\n');
    fprintf(fid,'Arquivo: %s\nFs: %.10g Hz\nAmostras: %d\nDuracao util: %.6f s (Tempo inicial descartado: %.1f s)\n',inputFile,Fs,nAmostras,tempo(end)-tempo(1),cfg.tInicioSeg);
    fprintf(fid,'Escala de amplitude: picoTesla (pT) [direta/sem conversao]\n');
    fprintf(fid,'Faixa espectral analisada: %.6g-%.6g cpm (%.8g-%.8g Hz)\n',cfg.cpmMin,cfg.cpmMax,cfg.freqMinHz,cfg.freqMaxHz);
    fprintf(fid,'Canais validos: %s\nCanais excluidos: %s\n\n',num2str(canaisValidosIdx),num2str(canaisExcluidosIdx));
    fprintf(fid,'RESULTADO GLOBAL PRINCIPAL\n');
    fprintf(fid,'Pico FFT global: %.8f cpm\n',fftInfo.picoCpm);
    fprintf(fid,'Frequencia dominante Welch: %.8f cpm (%.8f Hz)\n',welchInfo.picoCpm,welchInfo.picoCpm/60);
    fprintf(fid,'Classificacao global: %s\n',classeGlobal);
    fprintf(fid,'RMS global (pT): %.10g\nAmplitude pico a pico global (pT): %.10g\n',rmsGlobal,ampPPGlobal);
    fprintf(fid,'Amplitude PP media janelas %.1f s (pT): %.10g\nAmplitude PP mediana janelas (pT): %.10g\n',cfg.janelaAmplitudeSeg,ampPPMedia,ampPPMediana);
    fprintf(fid,'Potencia total %.3f-%.3f cpm (pT^2): %.10g\nPotencia no pico Welch (pT^2/Hz): %.10g\n',cfg.cpmMin,cfg.cpmMax,potenciaTotal,welchInfo.potenciaPico);
    fprintf(fid,'FWHM: %.8f cpm | limites %.8f a %.8f cpm\n',fwhm(1),fwhm(2),fwhm(3));
    fprintf(fid,'Entropia espectral normalizada: %.8f\n',entropia);
    fprintf(fid,'Instability Coefficient (IC): %.8f (Ciclos: %d, Dur. Media: %.2f s, DesvPad: %.2f s)\n\n', ...
        icGlobalInfo.IC, icGlobalInfo.nCiclos, icGlobalInfo.meanDur, icGlobalInfo.stdDur);
    
    fprintf(fid,'FRACAO DE POTENCIA POR FAIXA (%%)\n');
    for k=1:4, fprintf(fid,'%s: %.6f\n',nomesFaixas{k},fracaoPotencia(k)); end
    fprintf(fid,'\nPREVALENCIA TEMPORAL RSA (%%)\n');
    for k=1:4, fprintf(fid,'%s: %.6f\n',nomesFaixas{k},rsa.prevalencia(k)); end
    fprintf(fid,'\nCONTROLE INTERNO POR CANAL\nCanal\tValido_pct\tRMS(pT)\tFFT_cpm\tClasse_FFT\tWelch_cpm\tClasse_Welch\tDelta\tConcordancia\tEntropia\tIC\tStatus\n');
    for ch=1:nCanais
        fprintf(fid,'%d\t%.6f\t%.10g\t%.8f\t%s\t%.8f\t%s\t%.8f\t%s\t%.8f\t%.8f\t%s\n', ...
            controle(ch).canal,controle(ch).valido_pct,controle(ch).rms,controle(ch).fft_cpm,controle(ch).classe_fft, ...
            controle(ch).welch_cpm,controle(ch).classe_welch,controle(ch).delta,controle(ch).concordancia,controle(ch).entropia,controle(ch).ic,controle(ch).status);
    end
    fclose(fid);

    save(fullfile(outDir,'analise_global_completa.mat'),'tempo','Fs','canaisValidosIdx','canaisExcluidosIdx','sinalGlobal','sinalGlobalMediana', ...
        'fftInfo','welchInfo','classeGlobal','potenciaTotal','potFaixas','fracaoPotencia','fwhm','entropia','icGlobalInfo','ampPPGlobal','ampPPMedia','ampPPMediana','rmsGlobal','rsa','controle','cfg');

    fprintf('\n=============================================\nANÁLISE GLOBAL CONCLUÍDA\n=============================================\n');
    fprintf('Canais válidos: %s\n',num2str(canaisValidosIdx));
    fprintf('Frequência dominante Welch: %.4f cpm\n',welchInfo.picoCpm);
    fprintf('Classificação: %s\n',classeGlobal);
    fprintf('Resultados e arquivos .DAT salvos em:\n  %s\n',outDir);

    %% Plots Individuais Opcionais
    while true
        r=input('\nPlotar canal individual com sinal + FFT + Welch? [numero | 0=encerrar]: ','s');
        if isempty(r) || strcmp(strtrim(r),'0'), break; end
        ch=str2double(r);
        if ~isfinite(ch)||ch<1||ch>nCanais||~ismember(ch,canaisValidosIdx)
            fprintf('Canal inválido ou excluído.\n'); continue;
        end
        plot_canal_individual(ch,tempo,Xtrat(:,ch),Fs,cfg);
    end
end

%% ================= FUNÇÃO AUXILIAR DE LEITURA COMPATÍVEL =================
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
        error('Não foi possível ler uma matriz numérica válida a partir do arquivo fornecido.');
    end
end

%% ================= FUNÇÃO COM CRUZETA PARA SELEÇÃO DE CORTE =================
function [tCortado, fig] = plot_sinal_bruto_previa_com_cruzeta(t, X, tSugerido)
    nCanais = size(X, 2);
    fig = figure('Color', 'w', 'Name', 'SELEÇÃO DE CORTE TEMPORAL COM CRUZETA', 'Units', 'normalized', 'OuterPosition', [0.05 0.05 0.9 0.9]);
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
    xlabel('Tempo (segundos)');
    ylabel('Canais (pT)');
    title({'INSPEÇÃO PRÉVIA - Todos os canais (Sinal Completo)', ...
           'CLIQUE COM A CRUZETA no ponto do gráfico para definir o tempo de início (ou pressione ENTER para usar 0s)'});
    grid on;
    xlim([t(1) t(end)]);
    drawnow;

    try
        [xClick, ~, button] = ginput(1);
        if ~isempty(xClick) && ~isempty(button) && button == 1
            tCortado = max(0, xClick);
            line([tCortado tCortado], yl, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            text(tCortado + 2, yl(2)*0.9, sprintf('Corte Selecionado: %.2f s', tCortado), ...
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

%% ================= FUNÇÕES AUXILIARES PARA O IC =================
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

%% ================= OUTRAS FUNÇÕES AUXILIARES =================
function n=maior_sequencia_true(v)
    n=0; atual=0;
    for i=1:numel(v)
        if v(i), atual=atual+1; n=max(n,atual); else, atual=0; end
    end
end

function lista=parse_channel_list(txt,nMax)
    txt=strrep(txt,',',' '); txt=strrep(txt,'[',''); txt=strrep(txt,']','');
    if isempty(strtrim(txt)), lista=[]; return; end
    lista=unique(sscanf(txt,'%d').');
    if any(lista<1|lista>nMax), error('Há canais fora do intervalo permitido.'); end
end

function fig=plot_inspecao(t,X,fd,fc,fp,fr)
    n=size(X,2); fig=figure('Color','w','Name','INSPEÇÃO GLOBAL DE QUALIDADE','NumberTitle','off','Units','normalized','OuterPosition',[0 0 1 1]); hold on;
    bases=zeros(1,n); rot=cell(1,n);
    for ch=1:n
        x=X(:,ch); good=isfinite(x);
        if nnz(good)>=2
            if ~all(good), idx=(1:numel(x))'; x(~good)=interp1(idx(good),x(good),idx(~good),'linear','extrap'); end
            x=x-mean(x); esc=percentile_compat(x,95)-percentile_compat(x,5); if ~isfinite(esc)||esc<=eps, esc=std(x); end
            if isfinite(esc)&&esc>eps, x=x/esc; else, x(:)=NaN; end
        else, x(:)=NaN; end
        base=(n-ch)*2; bases(ch)=base; plot(t/60,x+base,'LineWidth',0.75);
        av={}; if fd(ch),av{end+1}='DADOS';end; if fc(ch),av{end+1}='CONSTANTE';end; if fp(ch),av{end+1}='PLATÔ';end; if fr(ch),av{end+1}='RMS BAIXO';end
        if isempty(av), rot{ch}=sprintf('Ch %d',ch); else, rot{ch}=sprintf('Ch %d [%s]',ch,strjoin_compatible(av,' | ')); end
    end
    [basesOrd, ordemTicks] = sort(bases, 'ascend');
    rotOrd = rot(ordemTicks);
    set(gca,'YTick',basesOrd,'YTickLabel',rotOrd); xlabel('Tempo (min)'); ylabel('Canais (normalizados apenas para inspeção)');
    title({'INSPEÇÃO DE QUALIDADE DOS CANAIS (PÓS-CORTE)','Feche esta janela para informar canais a excluir'}); grid on; xlim([t(1) t(end)]/60);
end

function fig=plot_empilhados(t,X,idx)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); hold on; n=numel(idx); amps=zeros(1,n);
    for k=1:n, x=X(:,idx(k)); amps(k)=percentile_compat(x,95)-percentile_compat(x,5); end
    esc=median(amps(isfinite(amps)&amps>0)); if isempty(esc), esc=max(std(X(:,idx),0,1)); end; if ~isfinite(esc)||esc<=0,esc=1;end
    off=1.6*esc; bases=zeros(1,n);
    for k=1:n, bases(k)=(n-k)*off; plot(t/60,X(:,idx(k))+bases(k),'LineWidth',0.75); end
    rotEmp = arrayfun(@(c)sprintf('Ch %d',c),idx,'UniformOutput',false);
    [basesOrd, ordemTicks] = sort(bases, 'ascend');
    rotEmpOrd = rotEmp(ordemTicks);
    set(gca,'YTick',basesOrd,'YTickLabel',rotEmpOrd); xlabel('Tempo (min)'); ylabel('Canais filtrados (pT)'); title(sprintf('Sinais pós-pipeline empilhados - %d canais usados',n)); grid on; xlim([t(1) t(end)]/60);
end

function info=calc_fft(x,Fs,cpmMin,cpmMax,nfftMin)
    x=x(:); nfft=max(nfftMin,2^nextpow2(numel(x))); w=hann_compat(numel(x)); Y=fft(x.*w,nfft); npos=floor(nfft/2)+1; Y=Y(1:npos); fHz=(0:npos-1)'*Fs/nfft;
    amp=abs(Y)/max(sum(w),eps); if numel(amp)>2, amp(2:end-1)=2*amp(2:end-1); end
    fCpm=fHz*60; m=fCpm>=cpmMin & fCpm<=cpmMax; fl=fCpm(m); al=amp(m); [~,ix]=max(al); pico=refinar_pico(fl,al,ix); ap=interp1(fl,al,pico,'linear');
    info=struct('fHz',fHz,'fCpm',fCpm,'amplitude',amp,'fLocalCpm',fl,'aLocal',al,'picoCpm',pico,'amplitudePico',ap);
end

function info=calc_welch(x,Fs,cfg)
    n=numel(x); wn=min(round(cfg.janelaWelchSeg*Fs),n); wn=max(min(wn,n),min(16,n)); w=hann_compat(wn); nov=min(wn-1,floor(cfg.sobreposicaoWelch*wn)); nfft=max(cfg.nfftWelch,2^nextpow2(wn));
    [pxx,fHz]=pwelch(x,w,nov,nfft,Fs,'psd'); fCpm=fHz*60; m=fCpm>=cfg.cpmMin & fCpm<=cfg.cpmMax; fl=fCpm(m); pl=pxx(m); [~,ix]=max(pl); pico=refinar_pico(fl,pl,ix); pp=interp1(fl,pl,pico,'linear');
    info=struct('fHz',fHz,'fCpm',fCpm,'pxx',pxx,'fLocalCpm',fl,'pLocal',pl,'idxPicoLocal',ix,'picoCpm',pico,'potenciaPico',pp,'janelaN',wn,'noverlap',nov,'nfft',nfft);
end

function pico=refinar_pico(f,y,ix)
    pico=f(ix); if ix>1 && ix<numel(y) && numel(f)>=2
        y1=log(max(y(ix-1),realmin)); y2=log(max(y(ix),realmin)); y3=log(max(y(ix+1),realmin)); den=y1-2*y2+y3;
        if isfinite(den)&&abs(den)>eps, d=0.5*(y1-y3)/den; d=max(-1,min(1,d)); pico=f(ix)+d*(f(2)-f(1)); end
    end
end

function val=integrar_banda(f,p,faixa)
    m=f>=faixa(1)&f<=faixa(2); if nnz(m)>=2,val=trapz(f(m),p(m));elseif nnz(m)==1,val=p(m);else,val=0;end
end

function f=calc_fwhm(freq,p,ix)
    meia=p(ix)/2; fe=NaN; fd=NaN;
    iesq=find(p(1:ix)<=meia,1,'last'); if ~isempty(iesq)&&iesq<ix, fe=interp_cross(freq(iesq),p(iesq),freq(iesq+1),p(iesq+1),meia); end
    idr0=find(p(ix:end)<=meia,1,'first'); if ~isempty(idr0), idr=ix+idr0-1; if idr>ix, fd=interp_cross(freq(idr-1),p(idr-1),freq(idr),p(idr),meia); end, end
    if isfinite(fe)&&isfinite(fd)&&fd>fe, largura=fd-fe; else, largura=NaN; end
    f=[largura fe fd];
end

function x=interp_cross(f1,p1,f2,p2,a), if abs(p2-p1)>eps,x=f1+(a-p1)*(f2-f1)/(p2-p1);else,x=(f1+f2)/2;end, end

function h=entropia_espectral(p)
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
        [~,ix]=max(pl); fd(k)=refinar_pico(fl,pl,ix); cl(k)=codigo_faixa(fd(k)); tm(k)=mean(t(ii))/60;
    end
    prev=nan(1,4); v=isfinite(cl); if any(v),for j=1:4,prev(j)=100*mean(cl(v)==j);end,end
    rsa=struct('tempoMin',tm,'freqCpm',freq,'potNorm',potNorm,'freqDomCpm',fd,'classe',cl,'prevalencia',prev,'janelaN',wn,'passoN',passo,'nfft',nfft);
end

function c=codigo_faixa(f)
    if f>=1&&f<2,c=1;elseif f>=2&&f<=4,c=2;elseif f>4&&f<9,c=3;elseif f>=9&&f<=12,c=4;else,c=NaN;end
end

function s=classificar_freq(f)
    c=codigo_faixa(f); nomes={'Bradigastria','Normogastria','Taquigastria','Banda intestinal'}; if isfinite(c),s=nomes{c};else,s='Fora das faixas fisiológicas';end
end

function fig=plot_fft_welch(fi,wi,cfg,classe,fwhm)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); subplot(2,1,1); plot(fi.fLocalCpm,fi.aLocal,'LineWidth',1.1); hold on; plot(fi.picoCpm,fi.amplitudePico,'o'); grid on; xlim([cfg.cpmMin cfg.cpmMax]); xlabel('Frequência (cpm)');ylabel('Amplitude (pT)');title(sprintf('FFT global - pico %.3f cpm',fi.picoCpm));
    subplot(2,1,2); plot(wi.fLocalCpm,wi.pLocal,'LineWidth',1.1);hold on;plot(wi.picoCpm,wi.potenciaPico,'o'); if all(isfinite(fwhm(2:3))), yl=ylim; plot([fwhm(2) fwhm(2)],yl,'--'); plot([fwhm(3) fwhm(3)],yl,'--');end; grid on;xlim([cfg.cpmMin cfg.cpmMax]);xlabel('Frequência (cpm)');ylabel('PSD (pT^2/Hz)');title(sprintf('PSD Welch global - pico %.3f cpm - %s',wi.picoCpm,classe));
end

function fig=plot_rsa(t,x,rsa,cfg)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); subplot(2,1,1); imagesc(rsa.tempoMin,rsa.freqCpm,rsa.potNorm); axis xy; colorbar; hold on; for y=[2 4 9], plot(xlim,[y y],'--'); end; ylim([cfg.cpmMin cfg.cpmMax]); ylabel('Frequência (cpm)'); title('RSA global normalizada por janela');
    subplot(2,1,2); plot(t/60,x,'LineWidth',0.8); grid on; xlabel('Tempo (min)');ylabel('Amplitude (pT)');title('Média global dos canais válidos');xlim([t(1) t(end)]/60);
end

function fig=plot_classificacao(frac,prev)
    fig=figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); nomes={'Bradi','Normo','Taqui','Intestinal'}; subplot(2,1,1);bar(frac);set(gca,'XTick',1:4,'XTickLabel',nomes);ylabel('Fração de potência (%)');title('Distribuição da potência global');grid on; subplot(2,1,2);bar(prev);set(gca,'XTick',1:4,'XTickLabel',nomes);ylabel('Prevalência RSA (%)');title('Prevalência temporal global');grid on;
end

function plot_canal_individual(ch,t,x,Fs,cfg)
    fi=calc_fft(x,Fs,cfg.cpmMin,cfg.cpmMax,cfg.nfftFFT); wi=calc_welch(x,Fs,cfg); cl=classificar_freq(wi.picoCpm); figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1]); subplot(3,1,1);plot(t/60,x);grid on;title(sprintf('Canal %d - sinal temporal',ch));xlabel('Tempo (min)');ylabel('Amplitude (pT)'); xlim([t(1) t(end)]/60); subplot(3,1,2);plot(fi.fLocalCpm,fi.aLocal);hold on;plot(fi.picoCpm,fi.amplitudePico,'o');grid on;xlim([cfg.cpmMin cfg.cpmMax]);title(sprintf('FFT canal %d - pico %.3f cpm',ch,fi.picoCpm));xlabel('Frequência (cpm)');ylabel('Amplitude (pT)'); subplot(3,1,3);plot(wi.fLocalCpm,wi.pLocal);hold on;plot(wi.picoCpm,wi.potenciaPico,'o');grid on;xlim([cfg.cpmMin cfg.cpmMax]);title(sprintf('Welch canal %d - pico %.3f cpm - %s',ch,wi.picoCpm,cl));xlabel('Frequência (cpm)');ylabel('PSD (pT^2/Hz)');
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