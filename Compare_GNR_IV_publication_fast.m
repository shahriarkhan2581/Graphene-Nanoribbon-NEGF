%% Compare_GNR_IV_publication_fast.m
% Fast/publication switch for transmission and ballistic I-V comparison.
% Required: H_7AGNR.m, H_11AGNR.m, H_8ZGNR.m, GNR_Current.m

clear; clc; close all;

%% Run mode
publicationMode = false; % false = fast test; true = final high-resolution run
useCache = true;         % reuse saved NEGF results when settings are unchanged

%% Physical settings
Ef   = 0;      % eV
Temp = 300;    % K
Vmax = 2.0;    % V
eta = 2e-3;    % eV; stable and sufficiently small for these plots

if publicationMode
    Nv = 101;
    NE = 1201;
    Erange = [-3 3];
    outTag = 'publication';
else
    Nv = 51;
    NE = 301;
    Erange = [-2.2 2.2];
    outTag = 'quick';
end

q = 1.602176634e-19;
h = 6.62607015e-34;
G0 = 2*q^2/h;

cacheFile = sprintf('GNR_NEGF_cache_%s.mat',outTag);

%% Calculate or load cached results
if useCache && isfile(cacheFile)
    fprintf('Loading cached NEGF data from %s...\n',cacheFile);
    S = load(cacheFile);
    Vall = S.Vall; I7raw = S.I7raw; I11raw = S.I11raw; I8raw = S.I8raw;
    T7 = S.T7; T11 = S.T11; T8 = S.T8; E = S.E;
else
    fprintf('Running 7-AGNR...\n');
    [H00,H01] = H_7AGNR();
    [Vall,I7raw,T7,E] = GNR_Current(H00,H01, ...
        'Ef',Ef,'Vmax',Vmax,'Nv',Nv,'NE',NE, ...
        'Erange',Erange,'eta',eta,'T',Temp);

    fprintf('Running 11-AGNR...\n');
    [H00,H01] = H_11AGNR();
    [V11,I11raw,T11,E11] = GNR_Current(H00,H01, ...
        'Ef',Ef,'Vmax',Vmax,'Nv',Nv,'NE',NE, ...
        'Erange',Erange,'eta',eta,'T',Temp);

    fprintf('Running 8-ZGNR...\n');
    [H00,H01] = H_8ZGNR();
    [V8,I8raw,T8,E8] = GNR_Current(H00,H01, ...
        'Ef',Ef,'Vmax',Vmax,'Nv',Nv,'NE',NE, ...
        'Erange',Erange,'eta',eta,'T',Temp);

    assert(max(abs(Vall-V11)) < 1e-12 && max(abs(Vall-V8)) < 1e-12, ...
        'Bias grids are inconsistent.');
    assert(max(abs(E-E11)) < 1e-12 && max(abs(E-E8)) < 1e-12, ...
        'Energy grids are inconsistent.');

    if useCache
        save(cacheFile,'Vall','I7raw','I11raw','I8raw','T7','T11','T8','E');
        fprintf('Saved cache: %s\n',cacheFile);
    end
end

%% Convert and clean
I7 = G0*I7raw;
I11 = G0*I11raw;
I8 = G0*I8raw;
T7 = max(real(T7(:)),0);
T11 = max(real(T11(:)),0);
T8 = max(real(T8(:)),0);

pos = Vall >= -1e-12;
V = Vall(pos);
I7p = I7(pos);
I11p = I11(pos);
I8p = I8(pos);

%% Plot style
c7  = [0.0000 0.2860 0.5720];
c11 = [0.8500 0.1500 0.1200];
c8  = [0.1000 0.5500 0.2500];

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 18.2 7.8],'PaperPositionMode','auto');
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

%% Transmission
ax1 = nexttile(tl,1);
plot(ax1,E,T8,':','Color',c8,'LineWidth',2.1); hold(ax1,'on');
plot(ax1,E,T11,'--','Color',c11,'LineWidth',1.8);
plot(ax1,E,T7,'-','Color',c7,'LineWidth',1.8);
xline(ax1,0,':','Color',[0.35 0.35 0.35], ...
    'LineWidth',0.8,'HandleVisibility','off');
xlabel(ax1,'Energy, $E$ (eV)','Interpreter','latex');
ylabel(ax1,'Transmission, $T(E)$','Interpreter','latex');
title(ax1,'(a) Transmission spectrum','FontWeight','normal');
legend(ax1,{'8-ZGNR','11-AGNR','7-AGNR'}, ...
    'Location','northwest','Box','off','FontSize',8);
xlim(ax1,Erange);
yMaxT = max([T7;T11;T8]);
if yMaxT <= 0, yMaxT = 1; end
ylim(ax1,[0,max(1,ceil(1.03*yMaxT))]);
grid(ax1,'on'); box(ax1,'on');

%% I-V
ax2 = nexttile(tl,2);
plot(ax2,V,I8p*1e6,':','Color',c8,'LineWidth',2.2); hold(ax2,'on');
markerStep = max(1,round(numel(V)/15));
plot(ax2,V,I11p*1e6,'--o','Color',c11,'LineWidth',1.6, ...
    'MarkerSize',3.4,'MarkerIndices',1:markerStep:numel(V), ...
    'MarkerFaceColor','w');
plot(ax2,V,I7p*1e6,'-','Color',c7,'LineWidth',1.8);
xlabel(ax2,'Bias voltage, $V$ (V)','Interpreter','latex');
ylabel(ax2,'Current, $I$ ($\mu$A)','Interpreter','latex');
title(ax2,'(b) Ballistic $I$-$V$ characteristics', ...
    'Interpreter','latex','FontWeight','normal');
legend(ax2,{'8-ZGNR','11-AGNR','7-AGNR'}, ...
    'Location','northwest','Box','off','FontSize',8);
xlim(ax2,[0 Vmax]);
yMaxI = max([I7p;I11p;I8p])*1e6;
if yMaxI <= 0, yMaxI = 1; end
ylim(ax2,[0,1.05*yMaxI]);
grid(ax2,'on'); box(ax2,'on');

set([ax1 ax2],'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.8,'TickDir','in','TickLength',[0.018 0.018], ...
    'XMinorTick','on','YMinorTick','on','GridAlpha',0.16,'Layer','top');

%% Diagnostics
[~,iE0] = min(abs(E-Ef));
fprintf('\n===== Transmission at E_F = %.3f eV =====\n',Ef);
fprintf('7-AGNR  : %.6f\n',T7(iE0));
fprintf('11-AGNR : %.6f\n',T11(iE0));
fprintf('8-ZGNR  : %.6f\n',T8(iE0));

fprintf('\n===== Current at V = %.2f V =====\n',V(end));
fprintf('7-AGNR  : %.6f microA\n',I7p(end)*1e6);
fprintf('11-AGNR : %.6f microA\n',I11p(end)*1e6);
fprintf('8-ZGNR  : %.6f microA\n',I8p(end)*1e6);

%% Save data and output
writetable(table(V,I7p,I11p,I8p, ...
    'VariableNames',{'Bias_V','I_7AGNR_A','I_11AGNR_A','I_8ZGNR_A'}), ...
    sprintf('GNR_IV_results_%s.csv',outTag));
writetable(table(E,T7,T11,T8, ...
    'VariableNames',{'Energy_eV','T_7AGNR','T_11AGNR','T_8ZGNR'}), ...
    sprintf('GNR_Transmission_results_%s.csv',outTag));

if publicationMode
    exportgraphics(fig,'GNR_NEGF_publication.png','Resolution',600);
    exportgraphics(fig,'GNR_NEGF_publication.pdf','ContentType','vector');
    savefig(fig,'GNR_NEGF_publication.fig');
    fprintf('\nSaved publication PNG, PDF, FIG, and CSV files.\n');
else
    exportgraphics(fig,'GNR_NEGF_quick.png','Resolution',200);
    fprintf('\nSaved quick PNG and CSV files.\n');
end
