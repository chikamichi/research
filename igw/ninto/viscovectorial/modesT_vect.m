% ----------------------------------------------
% IGW: propagation visqueuse totale ou partielle
% ----------------------------------------------

% {{{ pre-processing

clear all;
close all;

% flags

has_viscosity = 1;

% atmospheric data from USSA76 model
load modele1

% modele1:
%   altitude (0 - 500km)
%   density
%   density derivative
%   log density derivative

res     = 1;        % resolution: 1 by default, may be 10, 100... (diverges!)
binf    = 1;
bsup    = 9999;
Z       = modele1(binf:res:bsup,1);
rho     = modele1(binf:res:bsup,2);
drho    = modele1(binf:res:bsup,3);
dlgrho  = modele1(binf:res:bsup,4);

n_mod      = size(Z,1);
n_mod3     = 3*n_mod;   % for the sparse matrix
I(1:n_mod) = 1;         % plot helper

for j = 1:n_mod-1
  dz(j) = (Z(j+1)-Z(j))*1.e3;   % space step
end
dz(n_mod) = dz(n_mod-1);
dz = dz';

% velocity derivative
u0b = 0;      % bottom
u0t = 0;      % top
% linear gradient
alpha = (u0t-u0b)/((Z(n_mod)-Z(1))*1.e3);
alpha_local = alpha.*ones(n_mod,1);
u0 = Z.*(1.e3*alpha) + u0b;

% {{{ constants

M = 5.9768e+24;
G = 6.67e-11;
R = 6378.e3;

% constants }}}

% {{{ modal pre-processing

T  = [50 100 150].*60;
wn = 2*pi./T;

h  = 2500;
g  = 9.8;
c  = sqrt(g*h);
kx = wn./c;
ky = kx;

% modal pre-processing }}}

% {{{ iterators

nlat = max(size(ky));
nlon = max(size(kx));
nt   = max(size(wn));

% iterators }}}

if has_viscosity
    visc = 1.3e-5;  % TODO get this from data sets
end

load BestView2

% pre-processing }}}

% {{{ loop over k
for ik = 1:nlon % (nlon-1)/2 + 2:nlon
    ik;
    k = kx(ik);     % nombre d'onde horizontal courant

    %drawnow;        % DRAWNOW causes figure windows and their children to update and
                    %% flushes the system event queue. Any callbacks generated by incoming
                    %% events - e.g. mouse or key events - will be dispatched before
                    %% DRAWNOW returns.

    %figure(2*ik);   % FIGURE(H) makes H the current figure,  forces it to become visible,
                    %% and raises it above all other figures on the screen.  If Figure H
                    %% does not exist,  and H is an integer,  a new figure is created with
                    %% handle H.
    %set(2*ik, 'position', [44   250   500   750]);
    %set(gcf, 'Color', 'w');
    %set(gca, 'LineWidth', [2])
    %set(gca, 'fontsize', 18);
    %view(VV);
    %title(['Hwater = ', num2str(h), ' m and kx = ',  num2str(k), ' rad/m'], 'fontsize', 18);
    %%     xlabel('omega w (rad/s)', 'fontsize', 18);
    %xlabel('period T (min)', 'fontsize', 18);
    %ylabel('Vertical V_r_e_a_l (m/s)', 'fontsize', 18);
    %zlabel('altitude (km)', 'fontsize', 18);
    %grid on;
    %hold on;

    %figure(2*ik+1);
    %set(2*ik+1, 'position', [44   250   500   750]);
    %set(gcf, 'Color', 'w');
    %set(gca, 'LineWidth', [2])
    %set(gca, 'fontsize', 18);
    %view(VV);
    %title(['Hwater = ', num2str(h), ' m   and   kx = ',  num2str(k), ' rad/m'], 'fontsize', 18);
    %%    xlabel('omega (rad/s)', 'fontsize', 18);
    %xlabel('period T (min)', 'fontsize', 18);
    %ylabel('Vertical V_i_m_a_g (m/s)', 'fontsize', 18);
    %zlabel('altitude (km)', 'fontsize', 18);
    %grid on;
    %hold on;
    
    % loop over k }}}

    % {{{ boucle sur les trois frequences par mode, sachant qu'un mode est une altitude ici
    % avec nt = max(size(wn)) == 3 ici
    for iw = 1:nt %(nt-1)/2 + 2:nt
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vectorisation en K*X = B
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % {{{ constantes et utilitaires

        N2         = -g.*dlgrho;
        N2over2g   = (N2./(2*g)).^2;
        omega      = wn(iw);
        BigOmega   = omega - k.*u0;                % FIXME k a vectoriser quand je passerai en full 3D
        kxy2       = k^2;                          % FIXME no full 3D
        %phase      = k - omega                    % (pas unitaires) FIXME a modifier pour la full 3D
        phase      = k.*dz - 2*pi;
        expniphase = exp(-i.*phase);                % FIXME verifier cette formule
        %size(expniphase);
        kz2        = k.*k.*(N2./omega./omega - 1) - N2over2g;

        % pour le cas visqueux
        %BigGamma           = (BigOmega + i*visc*kxy2;
        % bis
        BigGamma           = BigOmega;
        BigGammadz         = BigGamma.*dz;
        BigGammadz2        = BigGammadz.*dz;
        phi                = 2.*visc.*sqrt(rho);
        phiOverBigGammadz2 = phi./BigGammadz2;

        %omega
        %kxy2
        %BigOmega(1)
        %BigOmega(end)
        %dz(end)
        %BigGamma(1)
        %BigGamma(end)
        %BigGammadz(1)
        %BigGammadz(end)
        %BigGammadz2(1)
        %BigGammadz2(end)
        %phi(1)
        %phi(end)
        %phiOverBigGammadz2(1)
        %phiOverBigGammadz2(end)
        %pause

        % constantes et utilitaires }}}

        % {{{ resolution vectorielle en differences finies

        % {{{ coefficients de propagation non visqueuse, hors des conditions initiales

        c1 = -(-k.*alpha.*2.*dz./BigOmega + dlgrho.*dz);
        c2 = -i.*k.*k.*2.*dz./BigOmega;
        c3 = -(-dlgrho.*dz);
        c4 = -(i.*2.*dz.*(BigOmega + dlgrho.*g./BigOmega)); % }}}
        
        % {{{ coefficients de propagation visqueuse, hors des conditions initiales
        
        % coefficients pour la perturbation de la vitesse horizontale
        %cu1 = 1 + i.*phiOverBigGammadz2;
        %cu2 =  i./2.*phiOverBigGammadz2;
        %cu3 = cu2;
        %cu4 = k./BigGamma;
        %for j = 3:n_mod-1
            %cu5(j) = (i*(u0(j+1) - u0(j-1)))/(2*BigGammadz(j));
        %end
        %% TODO FIXME a terme, il faudra peut-etre adopter cette technique du step-forward d'orde 1 pour les conditions au sommet pour d'autres coefficients ?
        %cu5(n_mod) = (i*(u0(n_mod) - u0(n_mod-1)))/(BigGammadz(n_mod));

        % bis
        cu1 =  1 + i.*phiOverBigGammadz2;
        cu2 = -i.*0.5.*phiOverBigGammadz2;
        cu3 =  cu2;
        cu4 = -k./BigGamma;
        for j = 3:n_mod-1
            cu5(j) = (i*(u0(j+1) - u0(j-1)))/(2*BigGammadz(j));
        end
        % TODO FIXME a terme, il faudra peut-etre adopter cette technique du step-forward d'orde 1 pour les conditions au sommet pour d'autres coefficients ?
        cu5(n_mod) = (i*(u0(n_mod) - u0(n_mod-1)))/(BigGammadz(n_mod));

        % coefficients pour la perturbation de la vitesse verticale
        %cw1 = -(dz.*dlgrho - 2.*k.*alpha./BigGamma);
        %cw2 =   2.*i.*kxy2.*dz./BigGamma;
        %cw3 = -(k.*phiOverBigGammadz2);
        %cw4 =   2.*k.*phiOverBigGammadz2;
        %cw5 = -(phiOverBigGammadz2);

        % bis
        cw1 = -(dz.*(dlgrho - 2.*k.*alpha./BigGamma));
        cw2 =  2.*i.*kxy2.*dz./BigGamma;
        cw3 = -(k.*phiOverBigGammadz2);
        cw4 = -2.*cw3;
        cw5 =  cw3;
        
        % coefficients pour la perturbation de pression
        %cp1 =  dz.*dlgrho;
        %cp2 =  2.*i.*(BigOmega + g./BigOmega.*dlgrho + i.*visc.*kxy2) + 2.*phi./dz;
        %cp3 = -phi./dz;
        %cp4 = -cp3; % }}}
        
        % bis
        cp1 =  dz.*dlgrho;
        cp2 =  2.*i.*dz.*(BigOmega + g./BigOmega.*dlgrho) + 2.*phi./dz;
        cp3 = -phi./dz;
        cp4 =  cp3; % }}}

        %isreal(cu1)
        %isreal(cu2)
        %isreal(cu3)
        %isreal(cu4)
        %isreal(cu5)
        %disp '---'
        %isreal(cw1)
        %isreal(cw2)
        %isreal(cw3)
        %isreal(cw4)
        %isreal(cw5)
        %disp '---'
        %isreal(cp1)
        %isreal(cp2)
        %isreal(cp3)
        %isreal(cp4)
        %pause

        % {{{ conditions initiales et coefficients associes

        % {{{ C.I. pures (1)

        % C.I.1 pour les perturbations de vitesse verticale et de pression
        % on se donne une perturbation unitaire, dans l'optique de l'etude des modes
        
        w(1)   = complex(1,0);

        kz2  = k*k*(-g/omega/omega*dlgrho(1) - 1) - 0.25/6.4e7;
        N2   = -g*dlgrho(1);
        if N2(1) >= omega*omega
            % cas d'une onde propagative
            p(1) = i/(k*k)*(alpha_local(1)*k - i*sqrt(kz2(1))*(omega - k*u0(1)) - 0.5*(omega - k*u0(1))*dlgrho(1))*w(1);
        else
            % cas d'une onde non propagative
            p(1) = i/(k*k)*(alpha_local(1)*k -   sqrt(kz2(1))*(omega - k*u0(1)) - 0.5*(omega - k*u0(1))*dlgrho(1))*w(1);
        end

        % C.I.1 pour la perturbation de la vitesse horizontale
        % FIXME je trouve un facteur 1/rho0 pour u(j)
        % u(j)   = 1/BigOmega(j) * (k*p(j) - i*alpha*w(j));
        % d'ou les coefficients :
        ui1  = i*alpha_local(1);
        ui2  = -k/BigOmega(1);
        % par ailleurs, ces valeurs explicites sont necessaires :
        % C.I.1 perturbation horizontale

        % C.I. pures (1) }}}

        % {{{ C.I. au premier pas de temps (2)

        b1  = -(-i*k*dz(1) * (1 / BigOmega(1) * (-i*alpha)) + 0.5*dlgrho(1)*dz(1) + 1);
        b2  = -(-i*k*dz(1) * (1 / BigOmega(1) * (k)));
        b3  = -(-0.5*dlgrho(1)*dz(1) + 1);
        b4  = -(i*BigOmega(1)*dz(1) + i/BigOmega(1)*dlgrho(1)*dz(1)*g);
        
        ui3 =  i*alpha_local(2)/BigOmega(2);
        ui4 = -k/BigOmega(2);

        name1      = 'runb.log';

        save_data1 = [mantexpnt(b1), mantexpnt(real(b2)), mantexpnt(imag(b2)), mantexpnt(b3), mantexpnt(real(b4)), mantexpnt(imag(b4))];

        %disp 'test!'
        %mantexpnt(min(real(cw3)))
        %mantexpnt(max(real(cw3)))
        %mantexpnt(real(cw3))
        %pause

        name2      = 'runc.log';
        save_data2 = [mantexpnt(c1), mantexpnt(real(c2)), mantexpnt(imag(c2)), mantexpnt(c3), mantexpnt(real(c4)), mantexpnt(imag(c4))];

        %size(real(cw1))
        %size(imag(cw1))
        %size(real(cw2))
        %size(imag(cw2))
        %size(real(cw3))
        %size(imag(cw3))
        %size(real(cw4))
        %size(imag(cw4))
        %size(real(cw5))
        %size(imag(cw5))

        name3      = 'runcw.log';
        save_data3 = [mantexpnt(cw1) mantexpnt(real(cw2)), mantexpnt(imag(cw2)), mantexpnt(real(cw3)), mantexpnt(imag(cw3)), mantexpnt(real(cw4)), mantexpnt(imag(cw4)), mantexpnt(real(cw5)), mantexpnt(imag(cw5))];
        %pause

        %size(real(cu1))
        %size(imag(cu1))
        %size(real(cu2))
        %size(imag(cu2))
        %size(real(cu3))
        %size(imag(cu3))
        %size(real(cu4))
        %size(imag(cu4))
        %size(real(cu5'))
        %size(imag(cu5'))

        name4      = 'runcu.log';
        save_data4 = [mantexpnt(real(cu1)), mantexpnt(imag(cu1)), mantexpnt(real(cu2)), mantexpnt(imag(cu2)), mantexpnt(real(cu3)), mantexpnt(imag(cu3)), mantexpnt(real(cu4)), mantexpnt(imag(cu4)), mantexpnt(cu5')];
        %pause

        %size(real(cp1))
        %size(imag(cp1))
        %size(real(cp2))
        %size(imag(cp2))
        %size(real(cp3))
        %size(imag(cp3))
        %size(real(cp4))
        %size(imag(cp4))

        name5      = 'runcp.log';
        save_data5 = [mantexpnt(cp1), mantexpnt(real(cp2)), mantexpnt(imag(cp2)), mantexpnt(cp3), mantexpnt(cp4)];
        %pause

        save(name1, 'save_data1', '-ascii');
        save(name2, 'save_data2', '-ascii');
        save(name3, 'save_data3', '-ascii');
        save(name4, 'save_data4', '-ascii');
        save(name5, 'save_data5', '-ascii');

        % C.I. au premier pas de temps (2) }}}

        % conditions initiales et coefficients associes }}}

        % {{{ construction de la matrice creuse et du second membre

        % construction des diagonales i.e. des vecteurs colonnes alternes, pour spdiags()

        % attention, il faut faire une permutation circulaire* pour les diagonales d'indices impairs, de facon a avoir le bon coefficient sur la premiere ligne. Ici, permutations d'ordre 3, selon le numero de colonne

       % construction : deux lignes, puis lecture globale (colonne par colonne, donc en alternant les valeurs -- ex. si une ligne de 1 et une ligne de 0, donnera une suite de 1 et 0 alternes)

        % {{{ construction du second membre

        % on construit le vecteur B des seconds membres (SM), avec :
        % - une partie conditions initiales
        % - TODO eventuellement, une partie de propagation non visqueuse (SM = 0)
        %   ce qui necessite de modifier K !
        % - une partie de propagation visqueuse avec les SM appropries

        % {{{ calculs des seconds membres

        % il y a deux sets de conditions limites, on commence donc au rang 3 pour matlab
        % le final est traite a part
        for j = 3:n_mod-1
            sm_tmp(j) = (phi(j)/(BigGamma(j)*dz(j)))*expniphase(j)*(u0(j+1) - 2*u0(j) + u0(j-1));
        end
        sm_tmp(n_mod) = (phi(n_mod)/(BigGamma(n_mod)*dz(n_mod-1)))*expniphase(n_mod)*(u0(n_mod) - 2*u0(n_mod-1) + u0(n_mod-2));
        sm_tmp = sm_tmp(3:end);
        
        sm_u   = i.*sm_tmp./2;
        sm_w   = k.*sm_tmp;
        sm_p   = zeros(1, n_mod-2);

        % vecteur des valeurs alternees
        sm = [sm_u'; sm_w'; sm_p'];
        sm = sm(:);

        name6      = 'runsm.log';
        save_data6 = [real(mantexpnt(sm_u)), imag(mantexpnt(sm_u)), mantexpnt(real(sm_w)), mantexpnt(imag(sm_w)), real(mantexpnt(sm_p)), imag(mantexpnt(sm_p))];
        save(name6, 'save_data6', '-ascii');
        disp 'coeff saved'
        pause

        %disp('sm :')
        %size(sm)
        %n_mod3-6

        % calculs des seconds membres }}}

        B = [0; w(1); p(1); ... % conditions initiales
             0;    0;    0; ... % le premier pas est non visqueux
                 sm         ... % la propagation visqueuse implique des termes non nuls,
            ];                  % tous calcules a priori avec les donnees du probleme

        %B
        %pause
        %disp('B :')
        %size(B)

        % construction du second membre }}}

       % {{{ FIXME diagonales pour le cas non visqueux
        %C13 = [c1'; c3'];
        %C2  = [zeros(1, size(c2,1)); c2'];             % *
        %C4  = [c4'; zeros(1, size(c4,1))];             % *

        %diagC13 = C13(:);
        %diagC2  = C2(:);
        %% FIXME pb de signe ici ! en attendant, je mets un moins pour retrouver
        %% des parties imaginaires positives...
        %diagC4  = -C4(:);
        
        %% matrice des conditions limites
        %cl0 = [1 0; 0 1];
        %cl1 = [b1 b2 1 0; b4 b3 0 1];
        %%cl  = [cl0 zeros(size(cl0,1), n_mod-size(cl0,1)) ; cl1 zeros(size(cl1,1), n_mod-size(cl1,1))]
        %n_mod2 = 2*n_mod;
        %cl  = [cl0 zeros(2, n_mod2-2) ; cl1 zeros(2, n_mod2-4)];

        %% construction de la matrice creuse des coefficients du systeme lineaire pour
        %% l'iteration en differences finies
        %% partie utile de la matrice diagonale a laquelle sont ajoutees les 2*2
        %% conditions limites/initiales
        
        %plop = ones(n_mod2, 1);
        %preK = spdiags([-plop diagC4 diagC13 diagC2 plop], 0:4, n_mod2, n_mod2);
        %preK = preK(1:n_mod2 - 4, :);
        %K    = [cl; preK];
        %%Kf = full(K);
        %%Kf(1:8,1:8)
        %%diagC2(1)
        %%diagC2(2)
        %%diagC2(3)
        %%diagC2(4)
        %%Kf(n_mod2-7:end,n_mod2-7:end)
        %%size(Kf)
        % }}}
        
        % {{{ diagonales pour le cas visqueux

        % helpers

        dimvisc = size(cu1, 1);
        myones  = ones( 1, dimvisc);
        myzeros = zeros(1, dimvisc);

        % diagonales alternees (les * denotent celles affectees d'une permutation circulaire)

        diag1   = [cw5';  cp4'; myzeros];       % *
        diag1   = diag1(:);

        diag2   = [cu2'; -myones; -myones];
        diag2   = diag2(:);

        diag3   = [myzeros myzeros myzeros];    % *
        diag3   = diag3';

        diag4   = [cw4'; cp2'; myzeros];        % *
        diag4   = diag4(:);

        diag5   = [cu1'; cw1'; cp1'];
        diag5   = diag5(:);

        diag6   = [myzeros; cu5; cw2'];        % *
        diag6   = diag6(:);

        diag7   = [cw3'; cp3'; cu4'];           % *
        diag7   = diag7(:);

        diag8   = [cu3'; myones; myones];
        diag8   = diag8(:);

        % cas visqueux }}}
        
        % {{{ matrices des conditions initiales

        ci0 = [1 ui1 ui2; ...
               0   1   0; ...
               0   0   1];

        ci1 = [0   0   0  1  ui3  ui4; ...
               0  b1  b2  0    1    0; ...
               0  b4  b3  0    0    1];

        ci  = [ci0 zeros(3, n_mod3-3) ; ci1 zeros(3, n_mod3-6)];

        % }}}

        % {{{ construction effective de la matrice creuse

        % construction de la matrice creuse des coefficients du systeme lineaire pour
        % l'iteration en differences finies
        % partie utile de la matrice diagonale a laquelle sont ajoutees les 2*2
        % conditions limites/initiales
        
        preK = spdiags([diag1 diag2 diag3 diag4 diag5 diag6 diag7 diag8], ...
                       -1:6, n_mod3, n_mod3);

        %precond_expnt = [10^mantexpnt(diag1)   ...
                         %10^mantexpnt(diag2)   ...
                         %10^mantexpnt(diag3)   ...
                         %10^mantexpnt(diag4)   ...
                         %10^mantexpnt(diag5)   ...
                         %10^mantexpnt(diag6)   ...
                         %10^mantexpnt(diag7)   ...
                         %10^mantexpnt(diag8)];
        size(diag1,1)
        %pause
        %precond = spdiags(repmat(precond_expnt, size(diag1,1), 1), -1:6, n_mod3, n_mod3);

        %full(preK)

        % troncature de la base, remplacee par les conditions initiales (6 lignes)
        
        preK    = preK(1:n_mod3 - 6, :);
        K       = [ci; preK];
        
        %disp 'cholesky'
        %pause
        %K = real(K);
        %B = real(B);
        %H = cholinc(K, '0');
        %spy(H)
        %pause
        %Kstar = (H')\K/H;
        %condKstar = condest(Kstar)
        %pause
        %[U, flag, relres, iter, resvec] = pcg(K, B, 1e-12, 2000, H', H);
        %iter
        %pause

        %precond = precond(1:n_mod3 - 6, :);
        %precond = [1 1 1       zeros(1, n_mod3 - 3) ; ...
                   %0 1 0       zeros(1, n_mod3 - 3) ; ...
                   %0 0 1       zeros(1, n_mod3 - 3) ; ...
                   %0 0 0 1 1 1 zeros(1, n_mod3 - 6) ; ...
                   %0 1 1 0 1 0 zeros(1, n_mod3 - 6) ; ...
                   %0 1 1 0 0 1 zeros(1, n_mod3 - 6) ; ...
                   %precond]';

        %size(diag(K,-1))
        %size(diag(precond,-1))
        %%full(diag(K,-1))
        %mantexpnt(diag(K,-1))
        %mantexpnt(diag(precond',-1))
        %mantexpnt(diag(K,-1)*diag(precond,-1)')
        %pause

        %disp('K :')
        %size(K)
        %Kf = full(K)
        %pause
        %figure
        %spy(K)
        %Kf(1:8,1:8)
        %diagC2(1)
        %diagC2(2)
        %diagC2(3)
        %diagC2(4)
        %Kf(n_mod2-7:end,n_mod2-7:end)
        %size(Kf)
  
        % construction de la matrice creuse }}}
 
        % construction de la matrice creuse }}}

        % {{{ resolution du systeme lineaire K*X = B

        %%%%%%%%%%%%%%%%%%%%
        U = K\B; % beautiful

        %U = pcg(K, B, [], 100, precond);
        %%%%%%%%%%%%%%%%%%%%

        %size(U)
        %U(1:13)
        %U(n_mod3-12:end)
        %p1

        % resolution de K*X = B }}}

        % {{{ resultats et calculs finaux

        % recuperation des resultats en w et p
        u = U(1:3:end); % TODO multiplier par precond_expnt ou un truc du genre !
        w = U(2:3:end);
        p = U(3:3:end);
        
        % calcul des kz explicites
        % FIXME coefficient a remplacer par la formule analytique, non ?
        kz2 = (k*k).*((-g/omega/omega).*dlgrho-1) - 0.25/6.4e7;
        
        % calcul des rho1 explicites
        rho1 = (-i.*dlgrho.*w)./BigOmega;

        % resultats et calculs finaux }}}

        % {{{ post-processing pour le pas courant

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
        axis([-2 2 0 max(Z)])
        v = axis;
        axis([v(1) v(2) Zmin Zmax])
        xlabel('Vertical V (m/s)', 'fontsize', 12);
        ylabel('Altitude (km)', 'fontsize', 16);
        title(['omega =  ', num2str(omega), ' rad/s'], 'fontsize', 12);
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
        title(['kx = ',  num2str(k), ' rad/m'], 'fontsize', 12);

        figure(100);
        F = getframe(gcf);
        [X, Map]  =  frame2im(F);
        file  =  ['Earth_Mode_L' num2str(round(1/k*1.e-3)) '_T_' num2str(2*pi/omega/60) '.tif'];
        imwrite(X, file, 'tif', 'Compression', 'none');

        disp('PRESS ENTER TO CONTINUE...');
        pause
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

        pause
        close

        figure(2*ik);
        plot3(I*2*pi/omega/60, real(w), Z, 'k');
        plot3(I*2*pi/omega/60, imag(w), Z, 'r');

        figure(2*ik+1);
        plot3(I*2*pi/omega/60, imag(w), Z, 'k');

        clear w u p rho1 kz2;

        % post-processing pour le pas courant }}}

    end
end     % resolution vectorielle en differences finies }}}

% boucle sur les trois frequences }}}

% {{{ post-processing

%Nsol = sqrt(N2(1));
%Nmax = max(sqrt(N2));
%Nmin = min(sqrt(N2));

%njump = 0;
%njump2 = size(wmax);

%n1 = 1;
%n2 = 1;

%h = figure;
%set(h, 'position', [2 700 800 600]);
%set(gcf, 'Color', 'w');
%font_size   =  22;
%font_sizeN  =  18;

%flag_plot   =  input('-Plot omega-k or T-lambda ? (0 or 1) ');

%subplot(2, 2, 1);
%if flag_plot  ==  0
    %pcolor(wn(n1:nt), kx(n2:nlon), log10(wmax(:,1+njump:njump2(2))));
%else
    %pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(wmax(:,1+njump:njump2(2))));
%end
%shading flat;
%v = axis;
%hold on;
%plot([Nmax Nmax], [v(3) v(4)], 'k');
%plot([Nsol Nsol], [v(3) v(4)], 'k');
%plot([Nmin Nmin], [v(3) v(4)], 'k');
%colorbar;
%title('W_r_e_a_l', 'fontsize', font_size);
%if flag_plot  ==  0
    %ylabel('kx (rad/m)', 'fontsize', font_size);
%else
    %ylabel('lambda (km)', 'fontsize', font_size);
%end
%set(gca, 'LineWidth', [2])
%set(gca, 'fontsize', font_sizeN);

%subplot(2, 2, 2);
%if flag_plot  ==  0
    %pcolor(wn(n1:nt), kx(n2:nlon), log10(umax(:,1+njump:njump2(2))));
%else
    %pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(umax(:,1+njump:njump2(2))));
%end
%shading flat;
%v = axis;
%hold on;
%plot([Nmax Nmax], [v(3) v(4)], 'k');
%plot([Nsol Nsol], [v(3) v(4)], 'k');
%plot([Nmin Nmin], [v(3) v(4)], 'k');
%colorbar;
%title('U_r_e_a_l', 'fontsize', font_size);
%set(gca, 'LineWidth', [2])
%set(gca, 'fontsize', font_sizeN);

%subplot(2, 2, 3);
%if flag_plot  ==  0
    %pcolor(wn(n1:nt), kx(n2:nlon), log10(pmax(:,1+njump:njump2(2))));
%else
    %pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(pmax(:,1+njump:njump2(2))));
%end
%shading flat;
%v = axis;
%hold on;
%plot([Nmax Nmax], [v(3) v(4)], 'k');
%plot([Nsol Nsol], [v(3) v(4)], 'k');
%plot([Nmin Nmin], [v(3) v(4)], 'k');
%colorbar;
%if flag_plot  ==  0
    %ylabel('kx (rad/m)', 'fontsize', font_size);
    %xlabel('omega (rad/s)', 'fontsize', font_size);
%else
    %ylabel('lambda (km)', 'fontsize', font_size);
    %xlabel('period T (min)', 'fontsize', font_size);
%end
%title('P_r_e_a_l', 'fontsize', font_size);
%set(gca, 'LineWidth', [2])
%set(gca, 'fontsize', font_sizeN);

%subplot(2, 2, 4);
%if flag_plot  ==  0
    %pcolor(wn(n1:nt), kx(n2:nlon), log10(rhomax(:,1+njump:njump2(2))));
%else
    %pcolor(2*pi./wn(n1:nt)/60, 2*pi./kx(n2:nlon)*1.e-3, log10(rhomax(:,1+njump:njump2(2))));
%end
%shading flat;
%v = axis;
%hold on;
%plot([Nmax Nmax], [v(3) v(4)], 'k');
%plot([Nsol Nsol], [v(3) v(4)], 'k');
%plot([Nmin Nmin], [v(3) v(4)], 'k');
%colorbar;
%if flag_plot  ==  0
    %xlabel('omega (rad/s)', 'fontsize', font_size);
%else
    %xlabel('period T (min)', 'fontsize', font_size);
%end
%title('Rho_r_e_a_l', 'fontsize', font_size);
%set(gca, 'LineWidth', [2])
%set(gca, 'fontsize', font_sizeN);

% post-processing }}}

