% TEST GRAVITY WAVES %

% version de jd, vectorisée

% pour le moment, il n'y a pas distinction entre les deux axes horizontaux,
% c'est du 3D planaire disons

% TODO ajouter un flag pour savoir si on calcule avec des L choisis ou bien des L dérivés des T via la vitesse de phase sismique
% TODO éventuellement, fichier ou vecteur de conf pour les arguments ?

function modesT_vect(p_u0b, p_u0t, p_T, p_L, p_res, plotflags, exportflags)
% arguments:
% - plotflags is [0|1, 0|1, 0|1, 0|1]
% - exportflags is [0|1, 0|1]
% plotflags are 0 by default, export_data is 1 (computations and export only)
% usage:
% - to generate data but no pics, modesT_vect(0, -50, (30:10:120), (50:50:400), 10, [0 0 0 0], [1 0]) eg.

disp '~~~~~~~~~~~~~~~~'
disp 'in modesT_vect()'
disp '~~~~~~~~~~~~~~~~'

% --- parameters handling and flags

%optargin = size(varargin, 2);
%stdargin = nargin - optargin;

disp '-----------'
disp 'en entrée :'
p_u0b
p_u0t
p_T
p_L

disp '---------------------------------'
disp 'après traitement des paramètres :'

if nargin < 1 || isempty(p_u0b)
    p_u0b = 0;
end

if nargin < 2 || isempty(p_u0t)
    p_u0t = 0;
end

if nargin < 3 || isempty(p_T)
    p_T = [10, 15, 20, 25, 30, 45, 50].*60;
else
    p_T = p_T.*60;
end
p_T

do_auto_L = 0;
if nargin < 4 || isempty(p_L)
    p_L = [50 100 150 200 250 300 350 400].*1e3;
else
    if isa(p_L, 'char')
        if p_L == 'auto'
            do_auto_L = 1; 
        else
            error('wrong value for wavelength (auto or an array)');
        end
    else
        p_L = p_L.*1e3;
    end
end
p_L
%pause

if nargin < 5 || isempty(p_res)
    p_res = 1
end

if nargin < 6 || isempty(plotflags)
    plotflags = [0, 0, 0, 0];
end

if nargin < 7 || isempty(exportflags)
    exportflags = [1, 0];
end

plot_cl          = plotflags(1);
plot_modes       = plotflags(2);
plot_amplitudes  = plotflags(3);
plot_scales      = plotflags(4);
do_export_data   = exportflags(1);
do_export_append = exportflags(2);

% modele1:
%   altitude (0 - 500km)
%   density
%   dérivée de la densité
%   dérivée du logarithme de la densité
load modele1

% tests à résolution variable
res     = p_res;    % resolution
binf    = 1;
bsup    = 9999;

Z       = modele1(binf:res:bsup,1);
rho     = modele1(binf:res:bsup,2);
drho    = modele1(binf:res:bsup,3);
dlgrho  = modele1(binf:res:bsup,4);

n_mod      = size(Z,1)
I(1:n_mod) = 1;

% calcul de la dérivée spatiale verticale
% TODO à vectoriser ^^ cf. version 3D avec les vents
for j = 1:n_mod-1
  dz(j) = (Z(j+1)-Z(j))*1.e3;     % pas spatial en mètres, utilisé dans les
end                               % dérivées en dz j'imagine
dz(n_mod) = dz(n_mod-1);          % FIXME ajouté pour la version vectorisée, bof

dz = dz';

% conditions limites sur la vitesse horizontale : pour le calcul du gradient de vent
u0b = p_u0b;      % vitesse en bas (bottom)
u0t = p_u0t;      % vitesse en haut (top)
% coefficient d'un gradient linéaire des vents égal à d(u0)/dz
alpha = (u0t-u0b)/((Z(n_mod)-Z(1))*1.e3);
%alpha
% profil des vitesses moyennes
u0 = Z.*(1.e3*alpha) + u0b;

M = 5.9768e+24;     % masse de la Terre
G = 6.67e-11;       % constante gravitationnelle
R = 6378.e3;        % rayon de la Terre

% {{{ traitement modal

disp '---------------------------'
T  = p_T

h  = 2500;
g  = 9.8;

wn = 2*pi./T;
%round(1/ikx*1.e-3)

if do_auto_L
    disp 'auto L'
    c  = sqrt(g*h);
    kx = wn./c
    ky = kx;
    L = round(2.*pi./kx.*1e-3)
else
    disp 'manual L'
    % essentiellement pour les tests à vent nul
    %L  = [10 15 20 25 50 100 150 200 250 300 350 400]*1.e3;
    kx = 2*pi./p_L
    ky = kx;
    L  = p_L./1e3
end

% }}}

% {{{ itérateurs spatio-temporels

nlat = max(size(ky));
nlon = max(size(kx));
nt   = max(size(wn));

% }}}

%----------------------------------
% Calculs principaux
% L'objectif est de calculer les modes pour les IGW, ie. la valeur des paramètres
% des équations du système à chaque altitude (« mode ») pour laquelle est fixé
% un kz approprié, sous la contrainte des conditions limites. La condition limite
% bottom fait office de condition initiale, couplée à une condition sur w, une
% condition sur P (continuité à l'interface).
%
% On considère trois périodes de temps, donc trois fréquences wn, donc trois
% nombre d'ondes k à travers la vitesse de phase imposée comme condition au limite
% implicite par continuité à l'interface. Il y a donc trois cas étudiés (ik).
% Pour chaque cas, on analyse les trois fréquences (iw). La boucle sur les iw
% résout le propagateur, sous les conditions initiales/limites.
%----------------------------------

load BestView2

% petits contrôles…

%pause
%disp '+++'
%binf
%res
%bsup
%T
%L
%disp '+++'
%pause

%----------------------------------

% boucle sur les nombres d'onde (sur les longeurs d'onde), ici trois.
% une figure par mode ik
for ik = 1:nlon % (nlon-1)/2 + 2:nlon
    
    %ik
    k = kx(ik);     % nombre d'onde horizontal courant
    l = L(ik);

    if plot_cl

        drawnow;        % DRAWNOW causes figure windows and their children to update and
        % flushes the system event queue. Any callbacks generated by incoming
        % events - e.g. mouse or key events - will be dispatched before
        % DRAWNOW returns.


        figure(2*ik);   % FIGURE(H) makes H the current figure,  forces it to become visible,
        % and raises it above all other figures on the screen.  If Figure H
        % does not exist,  and H is an integer,  a new figure is created with
        % handle H.

        set(2*ik, 'position', [44   250   500   750]);
        set(gcf, 'Color', 'w');
        set(gca, 'LineWidth', [2])
        set(gca, 'fontsize', 18);
        view(VV);
        title(['Hwater = ', num2str(h), 'm and kx = ',  num2str(k), ' rad/m'], 'fontsize', 18);
        %     xlabel('omega w (rad/s)', 'fontsize', 18);
        xlabel('period T (min)', 'fontsize', 18);
        ylabel('Vertical V_r_e_a_l (m/s)', 'fontsize', 18);
        zlabel('altitude (km)', 'fontsize', 18);
        grid on;
        hold on;

        figure(2*ik+1);
        set(2*ik+1, 'position', [44   250   500   750]);
        set(gcf, 'Color', 'w');
        set(gca, 'LineWidth', [2])
        set(gca, 'fontsize', 18);
        view(VV);
        title(['Hwater = ', num2str(h), 'm   and   kx = ',  num2str(k), ' rad/m'], 'fontsize', 18);
        %    xlabel('omega (rad/s)', 'fontsize', 18);
        xlabel('period T (min)', 'fontsize', 18);
        ylabel('Vertical V_i_m_a_g (m/s)', 'fontsize', 18);
        zlabel('altitude (km)', 'fontsize', 18);
        grid on;
        hold on;

    end

    % boucle sur les trois fréquences par mode, sachant qu'un mode est une altitude ici
    % avec nt = max(size(wn)) == 3 ici
    for iw = 1:nt %(nt-1)/2 + 2:nt
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vectorisation en K*U = 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        omega    = wn(iw);
        BigOmega = omega - k.*u0;

        % avec le propagateur en différences finies

        % calcul des coefficients pour les blocs diagonaux de la matrice creuse
        %c1 =  (i.*2.*dz.*k.*alpha)./BigOmega - dz.*dlgrho;
        %c2 =  (i*k*k)./BigOmega;
        %c3 =  dz.*dlgrho;
        %c4 = -(2*i).*dz.*(BigOmega + (g./BigOmega).*dlgrho);
        
        % tiré de la version de ninto, au signe près
        c1 = -(-k.*alpha.*2.*dz./BigOmega + dlgrho.*dz);
        c2 = -i.*k.*k.*2.*dz./BigOmega;
        c3 = -(-dlgrho.*dz);
        c4 = -(i.*2.*dz.*(BigOmega + dlgrho.*g./BigOmega));
        
        % condition initiale/limite sur P et w
        % calcul de la pression initiale sous contrainte de continuité
        w1  = complex(1,0);     % FIXME justification ?
        kz2 = k*k*(-g/omega/omega*dlgrho(1) - 1) - 0.25/6.4e7;
        N2  = -g*dlgrho(1);
        if N2 >= omega*omega
            % cas d'une onde propagative
            %disp 'onde propagative vent linéaire'
            p1 = i/(k*k)*(alpha*k - i*sqrt(kz2)*(omega - k*u0(1)) - 0.5*(omega - k*u0(1))*dlgrho(1))*w1;
        else
            % cas d'une onde non propagative
            p1 = i/(k*k)*(alpha*k -   sqrt(kz2)*(omega - k*u0(1)) - 0.5*(omega - k*u0(1))*dlgrho(1))*w1;
        end

        % calcul des coefficients pour les conditions initiales au step 1
        %b1 =  (i.*dz(1).*k.*alpha)./BigOmega(1) - 0.5.*dz(1).*dlgrho(1) - 1
        %b2 =  (i*k*k)./BigOmega(1)
        %b3 =  0.5.*dz(1).*dlgrho(1) - 1
        %b4 = -(2*i).*dz(1).*(BigOmega(1) + (g./BigOmega(1)).*dlgrho(1))
        
        % tiré de la version de ninto, au signe près
        b1 = -(-i*k*dz(1) * (1 / BigOmega(1) * (-i*alpha)) + 0.5*dlgrho(1)*dz(1) + 1);
        b2 = -(-i*k*dz(1) * (1 / BigOmega(1) * (k)));
        b3 = -(-0.5*dlgrho(1)*dz(1) + 1);
        b4 = -(i*BigOmega(1)*dz(1) + i/BigOmega(1)*dlgrho(1)*dz(1)*g);

        % vecteurs colonnes à 0 alternés pour spdiags()
        % attention, permutation circulaire pour les diagonales d'indices impairs*
        C13 = [c1'; c3'];
        C2  = [zeros(1, size(c2,1)); c2'];               % *
        C4  = [c4'; zeros(1, size(c4,1))];

        diagC13 = C13(:);
        %size(diagC13)
        diagC2  = C2(:);
        %size(diagC2)
        % FIXME pb de signe ici ! en attendant, je mets un moins pour retrouver
        % des parties imaginaires positives...
        diagC4  = -C4(:);
        %size(diagC4)
        
        % matrice des conditions limites
        cl0 = [1 0; 0 1];
        cl1 = [b1 b2 1 0; b4 b3 0 1];
        %cl  = [cl0 zeros(size(cl0,1), n_mod-size(cl0,1)) ; cl1 zeros(size(cl1,1), n_mod-size(cl1,1))]
        n_mod2 = 2*n_mod;
        cl  = [cl0 zeros(2, n_mod2-2) ; cl1 zeros(2, n_mod2-4)];

        % construction de la matrice creuse des coefficients du système linéaire pour
        % l'itération en différences finies
        % partie utile de la matrice diagonale à laquelle sont ajoutées les 2*2
        % conditions limites/initiales
        
        plop = ones(n_mod2, 1);
        preK = spdiags([-plop diagC4 diagC13 diagC2 plop], 0:4, n_mod2, n_mod2);
        preK = preK(1:n_mod2 - 4, :);
        %size(cl)
        %size(preK)
        K    = [cl; preK];
        %full(K)
        %pause

        % calcul du condition number, histoire de…
        %kn = svd(K)/svd(1/K);
        %disp('kn = ')
        %kn

        %Kf = full(K);
        %Kf(1:8,1:8)
        %diagC2(1)
        %diagC2(2)
        %diagC2(3)
        %diagC2(4)
        %Kf(n_mod2-7:end,n_mod2-7:end)
        %size(Kf)
        
        % résolution du système linéaire K*U = 0 -------------------------------
        % attention : faux pour la condition initiale qui n'est pas exprimée en
        % différences finies !
        % on construit le vecteur reste R nul partout sauf pour la CL
        R = [w1; p1; zeros(n_mod2-2, 1)];
        %pause
        U = K\R;
        %size(U)
        %U(1:10)
        %U(n_mod2-9:end)
        %p1

        % récupération des résultats en w et p
        w = U(1:2:end);
        %size(w)
        p = U(2:2:end);
        %size(p)
        
        % calcul du profil des vitesses intégrées explicites
        u = (-i.*alpha.*w + k.*p)./BigOmega;
        
        % calcul des kz explicites
        % FIXME coefficient à remplacer par la formule analytique, non ?
        kz2 = (k*k).*((-g/omega/omega).*dlgrho-1) - 0.25/6.4e7;
        
        % calcul des rho1 explicites
        rho1 = (-i.*dlgrho.*w)./BigOmega;

        % {{{ sauvegardes conditionnelles
        
        if do_export_data

            period = num2str(2*pi/omega/60);
            name   = ['data/wind/data_wind' num2str(u0b) 'to' num2str(u0t) '_L' num2str(l) '_T' period '.txt'];
            save_data = [real(w) imag(w) real(u) imag(u) real(p) imag(p) real(rho1) imag(rho1)];

            if do_export_append
                save(name, 'save_data', '-ascii', '-append');
            else
                save(name, 'save_data', '-ascii');
            end
        end
        
        % sauvegardes conditionnelles }}}
        
        wmax(ik, iw)    = max(abs(real(w)));
        umax(ik, iw)    = max(abs(real(u)));
        pmax(ik, iw)    = max(abs(real(p)));
        rhomax(ik, iw)  = max(abs(real(rho1)));
        KZ(ik, iw, :)   = sqrt(kz2);

        if max(abs(real(w))) >= 2.
            ctrw(ik, iw) = NaN;
        else
            ctrw(ik, iw) = 1.;
        end

        if sqrt(kz2(j)) >= 2*pi/50.e3
            A(ik, iw) = 1;
        else
            A(ik, iw) = NaN;
        end

        if plot_modes

            figure(100);
            set(100, 'position', [10 1 750 600])
            set(gcf, 'Color', 'w');
            fnt   =  20;
            fnt2  =  18;
            Zmax  =  max(Z);
            Zmin  =  min(Z);
            subplot(1, 4, 1);
            plot(real(w), Z, 'r', 'linewidth', 2);
            hold on;
            plot(imag(w), Z, 'g', 'linewidth', 2);
            hold off;
            %axis([-2 2 0 max(Z)])
            v = axis;
            axis([v(1) v(2) Zmin Zmax])
            xlabel('Vertical V (m/s)', 'fontsize', 12);
            ylabel('Altitude (km)', 'fontsize', 16);
            title(['Period = ', num2str(2*pi/omega/60), ' mn' ' (omega = ', num2str(omega), ' rad/s)'], 'fontsize', 12);
            subplot(1, 4, 2);
            plot(real(u), Z, 'r', 'linewidth', 2);
            hold on;
            plot(imag(u), Z, 'g', 'linewidth', 2);
            hold off;
            %axis([-20 20 0 max(Z)])
            v = axis;
            axis([v(1) v(2) Zmin Zmax]);
            xlabel('Horizontal V (m/s)', 'fontsize', 12);
            subplot(1, 4, 3);
            plot(real(p), Z, 'r', 'linewidth', 2);
            hold on;
            plot(imag(p), Z, 'g', 'linewidth', 2);
            hold off;
            %axis([-3000 3000 0 max(Z)])
            v = axis;
            axis([v(1) v(2) Zmin Zmax]);
            xlabel('Pressure (Pa)', 'fontsize', 12);
            subplot(1, 4, 4);
            plot(real(rho1), Z, 'r', 'linewidth', 2);
            hold on;
            plot(imag(rho1), Z, 'g', 'linewidth', 2);
            hold off;
            %axis([-0.1 0.1 0 max(Z)])
            v = axis;
            axis([v(1) v(2) Zmin Zmax]);
            xlabel('Density (kg/m3)', 'fontsize', 12);
            %title(['Lambda = ', num2str(round(1/k*1e-3)), ' km'  ' (kx = ',  num2str(k), ' rad/m)'], 'fontsize', 12);
            title(['Lambda = ', num2str(l/1e3), ' km'  ' (kx = ',  num2str(k), ' rad/m)'], 'fontsize', 12);

            figure(100);
            F = getframe(gcf);
            [X, Map]  =  frame2im(F);
            if u0t > 0
                wind = '+'
            else
                wind = ''
            end
            file  =  ['data/modes/mode_L' num2str(round(1/k*1.e-3)) '_T_' num2str(2*pi/omega/60) wind num2str(u0t) 'ms' '.tif'];
            imwrite(X, file, 'tif', 'Compression', 'none');

            %disp('PRESS ENTER TO CONTINUE...');
            %pause
            close

            figure(300);
            subplot(1, 3, 1);
            plot(real(sqrt(kz2)), Z, 'r');
            hold on;
            plot(imag(sqrt(kz2)), Z, 'g');
            xlabel('kz (rad/m)', 'fontsize', 12);
            ylabel('Altitude (km)', 'fontsize', 16);
            hold off;
            subplot(1, 3, 2);
            plot(sqrt(N2), Z, 'r');
            hold on;
            v = axis;
            plot([omega omega], [v(3) v(4)], 'g');
            hold off;
            xlabel('N and omega (rad/s)', 'fontsize', 12);
            subplot(1, 3, 3);
            plot(2*pi./real(sqrt(kz2))*1.e-3, Z, 'r');
            hold on;
            plot(2*pi./imag(sqrt(kz2))*1.e-3, Z, 'g');
            hold off;
            xlabel('lambda Z (km)', 'fontsize', 12);

            %pause
            close
        
        end

        if plot_cl

            figure(2*ik);
            plot3(I*2*pi/omega/60, real(w), Z, 'k');
            plot3(I*2*pi/omega/60, imag(w), Z, 'r');

            figure(2*ik+1);
            plot3(I*2*pi/omega/60, imag(w), Z, 'k');
        
        end

        clear w u p rho1 kz2
    end
end

wmax
%size(wmax)
%pause
%umax
%pause
%pmax
%pause
%rhomax
%pause

name1      = 'runb.log';
save_data1 = [b1, real(b2), imag(b2), b3, real(b4), imag(b4)];
name2      = 'runc.log';
save_data2 = [c1, real(c2), imag(c2), c3, real(c4), imag(c4)];

save(name1, 'save_data1', '-ascii');
save(name2, 'save_data2', '-ascii');

if plot_scales

    %Nsol = sqrt(N2(1));
    %Nmax = max(sqrt(N2));
    %Nmin = min(sqrt(N2));

    njump  = 0;
    njump2 = size(wmax);

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

    subplot(2, 2, 1);
    if flag_plot  ==  0
        pcolor(wn(n1:nt), kx(n2:nlon), log10(wmax(:,1+njump:njump2(2))));
    else
        size(wn)
        size(kx)
        size(wmax)
        pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(wmax(:,1+njump:njump2(2))));
        % wmax : toutes les lignes, ie. tous les ik (L) pour un certain nombre de iw (T), ici tous
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
        ylabel('kx (rad/m)', 'fontsize', font_size);
    else
        ylabel('lambda (km)', 'fontsize', font_size);
    end
    set(gca, 'LineWidth', [2])
    set(gca, 'fontsize', font_sizeN);

    subplot(2, 2, 2);
    if flag_plot  ==  0
        pcolor(wn(n1:nt), kx(n2:nlon), log10(umax(:,1+njump:njump2(2))));
    else
        pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(umax(:,1+njump:njump2(2))));
    end
    shading flat;
    %v = axis;
    %hold on;
    %plot([Nmax Nmax], [v(3) v(4)], 'k');
    %plot([Nsol Nsol], [v(3) v(4)], 'k');
    %plot([Nmin Nmin], [v(3) v(4)], 'k');
    colorbar;
    title('U_r_e_a_l', 'fontsize', font_size);
    set(gca, 'LineWidth', [2])
    set(gca, 'fontsize', font_sizeN);

    subplot(2, 2, 3);
    if flag_plot  ==  0
        pcolor(wn(n1:nt), kx(n2:nlon), log10(pmax(:,1+njump:njump2(2))));
    else
        pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(pmax(:,1+njump:njump2(2))));
    end
    shading flat;
    %v = axis;
    %hold on;
    %plot([Nmax Nmax], [v(3) v(4)], 'k');
    %plot([Nsol Nsol], [v(3) v(4)], 'k');
    %plot([Nmin Nmin], [v(3) v(4)], 'k');
    colorbar;
    if flag_plot  ==  0
        ylabel('kx (rad/m)', 'fontsize', font_size);
        xlabel('omega (rad/s)', 'fontsize', font_size);
    else
        ylabel('lambda (km)', 'fontsize', font_size);
        xlabel('period T (min)', 'fontsize', font_size);
    end
    title('P_r_e_a_l', 'fontsize', font_size);
    set(gca, 'LineWidth', [2])
    set(gca, 'fontsize', font_sizeN);

    subplot(2, 2, 4);
    if flag_plot  ==  0
        pcolor(wn(n1:nt), kx(n2:nlon), log10(rhomax(:,1+njump:njump2(2))));
    else
        pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(rhomax(:,1+njump:njump2(2))));
    end
    shading flat;
    %v = axis;
    %hold on;
    %plot([Nmax Nmax], [v(3) v(4)], 'k');
    %plot([Nsol Nsol], [v(3) v(4)], 'k');
    %plot([Nmin Nmin], [v(3) v(4)], 'k');
    colorbar;
    if flag_plot  ==  0
        xlabel('omega (rad/s)', 'fontsize', font_size);
    else
        xlabel('period T (min)', 'fontsize', font_size);
    end
    title('Rho_r_e_a_l', 'fontsize', font_size);
    set(gca, 'LineWidth', [2])
    set(gca, 'fontsize', font_sizeN);

end

