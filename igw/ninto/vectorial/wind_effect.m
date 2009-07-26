function wind_effect(u0b, u0t, T, L)

%clear all;
%close all;

% TODO s'assurer de l'existence des dossiers, les crééer sinon
% TODO indiquer la res dans les noms de fichiers, ça éviterait de les ouvrir…
%      trouver comment splitter la chaîne
% TODO bon ça ne gère pas les galères liées aux bornes, qui faussent la res
% TODO virer les output graphiques à chaque nouveau run

% {{{ flags

global force_compute

force_compute   = 0;
check_res       = 0;  % eg. assume res == 1
plot_amplitudes = 0;
plot_pcolors    = 1;

% }}}

% {{{ parameters handling

if u0b == 0 && u0t == 0
    error('at least one of the wind value must be different from zero');
end

if nargin < 1 || isempty(u0b)
    u0b = 0;
end

if nargin < 2 || isempty(u0t)
    u0t = 0;
end

if nargin < 3 || isempty(T)
    T = [15, 30, 60, 120];
end

if nargin < 4 || isempty(L)
    L  = [20 50 100 200 300 400];
else
    if L == 'auto'
        h  = 2500;
        g  = 9.8;
        c  = sqrt(g*h);
        wn = 2*pi./(T.*60)
        k  = wn./c
        L  = round(2.*pi./k.*1e-3)
        %pause
    else
        L = L.*1e3;
    end
end


% }}}

% {{{ import no-wind values sets
% these are obtained by running modesT_vect setting winds to 0 m/s

prefix = 'data/';
res  = 1;     % attention ! ne pas faire de tests à résolution inf. à 10, le code devient nécessairement instable
binf = 1;
bsup = 9999;

% USSA model
load modele1;
Z  = modele1(binf:res:bsup,1);
Zmax  =  max(Z);
Zmin  =  min(Z);

omega = 2.*pi./T;

if force_compute
    disp '--------------------------------------------------------------'
    disp 'force_compute is enabled: all data sets will be computed again'
    disp '--------------------------------------------------------------'
else
    disp '----------------------------------------------------------------------------------'
    disp 'force_compute is disabled: only outdated or unavailable data sets will be computed'
    disp '----------------------------------------------------------------------------------'
end
    
% {{{ getting refs

wind_ref = [];
do_L     = [];
do_T     = [];
do_check = [];

if force_compute
    do_L = L
    do_T = T
else
    % TODO FIXME ça ne gère que w pour l'instant, il faut des do_L et do_T à quatre composantes
    % TODO disp à remplacer par un compteur, n'afficher que les recomputés
    for j = 1:length(L)
        for k = 1:length(T)
            do_it    = 0;
            name = [prefix 'wind/data_wind0to0_L' num2str(L(j)) '_T' num2str(T(k)) '.txt'];
            if exist(name, 'file')
                fprintf('%s %s %s %s', '"',  name, '"', ' exists');
                if check_res
                    fprintf('%s', ', checking its resolution…');
                    do_check = load(name);
                    % force computation if the provided file has wrong resolution
                    %taille = size(do_check(:,1),1)
                    %bsup-binf/res+1
                    %pause
                    if size(do_check(:,1),1) ~= (bsup-binf/res+1)
                        fprintf('%s', ' must be computed again');
                        do_it = 1;
                    else
                        fprintf('%s', ' ok');
                    end
                    clear do_check;
                end
                fprintf('\n');
            else
                disp(['"' name '"' ' does not exists and will be computed'])
                do_it = 1;
            end

            if do_it
                do_L = [do_L L(j)];
                do_T = [do_T T(k)];
            end
        end
    end
    do_L = unique(do_L)
    do_T = unique(do_T)
end

if ~isempty(do_T) || ~isempty(do_L)
    disp '+++++++++'
    do_T
    do_L
    res
    disp '+++++++++'
    %pause
    modesT_vect(0, 0, do_T, do_L, res);  % outputs in data_wind0to0_L[L]_T[T].txt
end

for j = 1:length(L)
    for k = 1:length(T)
        name = [prefix 'wind/data_wind0to0_L' num2str(L(j)) '_T' num2str(T(k)) '.txt'];
        disp(['loading ' name]);
        wind_ref(:,:,j,k) = data_processing(name, 'real', binf, 1, bsup);
    end
end

%pause

% getting refs }}}

% check
%for j = 1:length(L)
    %for k = 1:length(T)
        %figure(k);
        %set(k, 'position', [10 1 750 600]);
        %set(gcf, 'Color', 'w');

        %plot(wind_ref(:,1,j,k), Z, 'r', 'linewidth', 2);
        %v = axis;
        %axis([v(1) v(2) Zmin Zmax]);
        %xlabel('Vertical V (m/s)', 'fontsize', 12);
        %ylabel('Altitude (km)', 'fontsize', 16);
        %title(['Period = ', num2str(T), ' mn' ' (omega = ', num2str(omega), ' rad/s)'], 'fontsize', 12);
        %pause
    %end
%end

% }}}

% {{{ computations for the input parameters

% {{{ compute only required data files

if force_compute
    do_L = L
    do_T = T
else
    % TODO FIXME ça ne gère que w pour l'instant, il faut des do_L et do_T à quatre composantes
    % TODO disp à remplacer par un compteur, n'afficher que les recomputés
    do_L     = [];
    do_T     = [];
    do_check = [];
    for j = 1:length(L)
        for k = 1:length(T)
            do_it = 0;
            name = [prefix 'wind/data_wind' num2str(u0b) 'to' num2str(u0t) '_L' num2str(L(j)) '_T' num2str(T(k)) '.txt'];
            if exist(name, 'file')
                fprintf('%s %s %s %s', '"',  name, '"', ' exists');
                if check_res
                    fprintf('%s', ', checking its resolution…');
                    do_check = load(name);
                    % force computation if the provided file has wrong resolution
                    if size(do_check(:,1),1) ~= (bsup-binf/res+1)
                        fprintf('%s', ' must be computed again');
                        do_it = 1;
                    else
                        fprintf('%s', ' ok');
                    end
                    clear do_check;
                end
                fprintf('\n');
            else
                disp(['"' name '"' ' does not exists and will be computed'])
                do_it = 1;
            end

            if do_it
                do_L = [do_L L(j)];
                do_T = [do_T T(k)];
            end
        end
    end
    do_L = unique(do_L)
    do_T = unique(do_T)
end
%pause

% compute only required data files }}}

if ~isempty(do_T) || ~isempty(do_L)
    disp '+++++++++'
    u0b
    u0t
    do_T
    do_L
    res
    disp '+++++++++'
    %pause
    modesT_vect(u0b, u0t, do_T, do_L, res);  % outputs in data_wind[u0b]to[u0t]_L[L]_T[T].txt
end

% computations }}}

% {{{ getting data

data = [];
for j = 1:length(L)
    for k = 1:length(T)
        name = [prefix 'wind/data_wind' num2str(u0b) 'to' num2str(u0t) '_L' num2str(L(j)) '_T' num2str(T(k)) '.txt'];
        % TODO à remplacer par un compteur
        disp(['loading ' name]);
        data(:,:,j,k) = data_processing(name, 'real', binf, 1, bsup);
    end
end

% getting data }}}

% {{{ amplitudes processing

% construction d'un vecteur des max abs des amplitudes,
% dans l'ordre des longueurs d'onde puis des périodes
% TODO FIXME il faut choper les squeeze de w pour tous les L et tous les T
% mais il faut réfléchir à comment ploter ça…
% idées :
% faire les deux, ie. pour chaque L, le max en fonction de T, puis pour chaque T, le max en fonction de L,
% et essayer un pcolor pour avoir le max pour tous les couples L/T (mais ça doit plus ou moins revenir au même que le pcolor du script principal : ses pcolor plottent l'amplification relative par rapport au sol, ici ce serait, intégré, par rapport à un cas de référence à vent nul, donc c'est pas exactement la même chose quand même)
% ne pas oublier de plotter le cas de référence à vent nul, d'ailleurs
% un truc qui serait chanmé : un pcolor en couche, avec une croix indiquant la position du cas de référence, une autre indiquant la position du cas étudié là, le vecteur vent moyen, en gros. À réfléchir…

% TODO boucler de 1 à 4 ligne 108 pour avoir w, u, p et r
% TODO à faire, pareil, pour la référence vent nul ! il faut un amp_max_ref

% a est data, b est amp_max, test est test ^^
% exemple de data :
%data(:,:,1,1) = [1 2; 3 4; 5 6];
%data(:,:,1,2) = [11 22; 33 44; 55 66];
%data(:,:,2,1) = [101 202; 303 404; 505 606];
%data(:,:,2,2) = [-1 -2; -3 -4; -5 -6];

disp 'processing des amplitudes'

amp_max_builder = [];
amp_ref_builder = [];
amp_max         = [];
amp_ref         = [];
amp_tmp         = squeeze(data(:,1,:,:));   % w real values 
amp_tmp2        = squeeze(wind_ref(:,1,:,:));
% attention, ça lit d'abord les j, puis les k : (1,1), (2,1), (1,2), (2,2)

test  = squeeze(min(amp_tmp) +max(amp_tmp))  < 0;   % look for the sign
test2 = squeeze(min(amp_tmp2)+max(amp_tmp2)) < 0;
% then process each w record
for j = 1:size(test,1)
    for k = 1:size(test,2)
        if test(j,k) == 1
            amp_max_builder = [amp_max_builder min(amp_tmp(:,j,k))];
        else
            amp_max_builder = [amp_max_builder max(amp_tmp(:,j,k))];
        end

        if test2(j,k) == 1
            amp_ref_builder = [amp_ref_builder min(amp_tmp2(:,j,k))];
        else
            amp_ref_builder = [amp_ref_builder max(amp_tmp2(:,j,k))];
        end
    end
    amp_max         = [amp_max; amp_max_builder];   % L in columns (one max per line for a given T),
    amp_ref         = [amp_ref; amp_ref_builder];   % T in rows (one max per column for a given L)
    amp_max_builder = [];
    amp_ref_builder = [];
end

% amplitudes processing }}}

% {{{ plot 'em all

% le squeeze est cool car il donnera autant de colonnes que de T, autant de lignes que de L
% ça permet de plotter direct versus L et T, pour un T ou un L donné respectivement !

% TODO plot plot !
% FIXME pas vraiment sûr que ça fonctionne, les valeurs de amp_max sont cheloues 
% tout vérifier pas à pas, ie. quelles valeurs ça manipule par rapport aux fichiers .txt
% (faire des load manuels et vérifier les amplitudes max etc.)

if plot_amplitudes

for k = 1:length(T)
    figure(2*k);
    hold on;
    plot(L, amp_max(:,k), 'r', L, amp_ref(:,k), '--b');
    title(['maximum w amplitude vs wave length, for period ', num2str(T(k))]);
    legend_wind = ['wind: ' num2str(u0b) ' to ' num2str(u0t) 'm/s'];
    legend(legend_wind, 'no wind (ref)');
    hold off;

    F = getframe(gcf);
    [X, Map]  =  frame2im(F);
    file  =  ['data/output/amplitudes/by-period/amp_T' num2str(T(k)) '_wind' num2str(u0b) 'to' num2str(u0t) 'ms' '.tif'];
    imwrite(X, file, 'tif', 'Compression', 'none');
end

close all;

for j = 1:length(L)
    figure(2*j+1);
    hold on;
    plot(T, amp_max(j,:), 'r', T, amp_ref(j,:), '--b');
    title(['maximum w amplitude vs period, for wave length ', num2str(L(j))]);
    legend_wind = ['wind: ' num2str(u0b) ' to ' num2str(u0t) 'm/s'];
    legend(legend_wind, 'no wind (ref)');
    hold off;

    F = getframe(gcf);
    [X, Map]  =  frame2im(F);
    file  =  ['data/output/amplitudes/by-wavelength/amp_L' num2str(L(j)) '_wind' num2str(u0b) 'to' num2str(u0t) 'ms' '.tif'];
    imwrite(X, file, 'tif', 'Compression', 'none');
end

end

if plot_pcolors

    %njump  = 0;
    %njump2 = size(amp_max);

    n1 = 1;
    n2 = 1;

    h = figure;
    set(h, 'position', [2 700 800 600]);
    set(gcf, 'Color', 'w');
    font_size  =  22;
    font_sizeN =  18;

    %flag_plot   =  input('-Plot omega-k or T-lambda ? (0 or 1) ');
    flag_plot = 1;
    colormap(hot);

    %subplot(2, 2, 1);
    if flag_plot  ==  0
        pcolor(wn, kx, amp_max-amp_ref);
    else
        %size(wn)
        %size(kx)
        %size(amp_max)
        %pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(wmax(:,1+njump:njump2(2))));
        % wmax : toutes les lignes, ie. tous les ik (L) pour un certain nombre de iw (T), ici tous
        pcolor(T, L, amp_max-amp_ref);
    end
    shading flat;
    %v = axis;
    %hold on;
    %plot([Nmax Nmax], [v(3) v(4)], 'k');
    %plot([Nsol Nsol], [v(3) v(4)], 'k');
    %plot([Nmin Nmin], [v(3) v(4)], 'k');
    colorbar;
    title('W_r_e_a_l', 'fontsize', font_size);
    if flag_plot  ==  0
        xlabel('omega (rad/s)', 'fontsize', font_size);
        ylabel('kx (rad/m)', 'fontsize', font_size);
    else
        xlabel('period T (min)', 'fontsize', font_size);
        ylabel('lambda (km)', 'fontsize', font_size);
    end
    set(gca, 'LineWidth', [2])
    set(gca, 'fontsize', font_sizeN);

end

%close all;

% plot 'em all }}}

% }}}

