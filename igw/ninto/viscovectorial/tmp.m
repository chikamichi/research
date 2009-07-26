% ----------------------------------------------
% IGW: propagation visqueuse totale ou partielle
% ----------------------------------------------

% {{{ pre-processing

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

    % figures...

    % loop over k }}}

    % {{{ loop over frequencies
    %     with nt = max(size(wn))
    for iw = 1:nt %(nt-1)/2 + 2:nt
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % solving K*U = B
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % {{{ constants, helpers

        N2         = -g.*dlgrho;
        N2over2g   = (N2./(2*g)).^2;
        omega      = wn(iw);
        BigOmega   = omega - k.*u0;
        kxy2       = k^2;
        phase      = k.*dz - 2*pi;
        expniphase = exp(-i.*phase);
        kz2        = k.*k.*(N2./omega./omega - 1) - N2over2g;

        BigGamma           = (BigOmega + i*visc*kxy2;
        BigGammadz         = BigGamma.*dz;
        BigGammadz2        = BigGammadz.*dz;
        phi                = 2.*visc.*sqrt(rho);
        phiOverBigGammadz2 = phi./BigGammadz2;

        % constants, helpers }}}

        % {{{ finite diff. solving

        % {{{ non viscuous

        c1 = -(-k.*alpha.*2.*dz./BigOmega + dlgrho.*dz);
        c2 = -i.*k.*k.*2.*dz./BigOmega;
        c3 = -(-dlgrho.*dz);
        c4 = -(i.*2.*dz.*(BigOmega + dlgrho.*g./BigOmega));
        
        % non viscuous }}}
        
        % {{{ viscuous
        
        coefficients pour la perturbation de la vitesse horizontale
        cu1 = 1 + i.*phiOverBigGammadz2;
        cu2 =  i./2.*phiOverBigGammadz2;
        cu3 = cu2;
        cu4 = k./BigGamma;
        for j = 3:n_mod-1
            cu5(j) = (i*(u0(j+1) - u0(j-1)))/(2*BigGammadz(j));
        end
        cu5(n_mod) = (i*(u0(n_mod) - u0(n_mod-1)))/(BigGammadz(n_mod));

        coefficients pour la perturbation de la vitesse verticale
        cw1 = -(dz.*dlgrho - 2.*k.*alpha./BigGamma);
        cw2 =   2.*i.*kxy2.*dz./BigGamma;
        cw3 = -(k.*phiOverBigGammadz2);
        cw4 =   2.*k.*phiOverBigGammadz2;
        cw5 = -(phiOverBigGammadz2);

        coefficients pour la perturbation de pression
        cp1 =  dz.*dlgrho;
        cp2 =  2.*i.*(BigOmega + g./BigOmega.*dlgrho + i.*visc.*kxy2) + 2.*phi./dz;
        cp3 = -phi./dz;
        cp4 = -cp3
        
        % viscuous }}}
        
        % finite diff. solving }}}

        % {{{ initial conditions

        % {{{ step 1

        w(1)   = complex(1,0);

        kz2  = k*k*(-g/omega/omega*dlgrho(1) - 1) - 0.25/6.4e7;
        N2   = -g*dlgrho(1);
        if N2(1) >= omega*omega
            p(1) = i/(k*k)*(alpha_local(1)*k - i*sqrt(kz2(1))*(omega - k*u0(1)) - 0.5*(omega - k*u0(1))*dlgrho(1))*w(1);
        else
            p(1) = i/(k*k)*(alpha_local(1)*k -   sqrt(kz2(1))*(omega - k*u0(1)) - 0.5*(omega - k*u0(1))*dlgrho(1))*w(1);
        end

        ui1  = i*alpha_local(1);
        ui2  = -k/BigOmega(1);

        % step 1 }}}

        % {{{ step 2

        b1  = -(-i*k*dz(1) * (1 / BigOmega(1) * (-i*alpha)) + 0.5*dlgrho(1)*dz(1) + 1);
        b2  = -(-i*k*dz(1) * (1 / BigOmega(1) * (k)));
        b3  = -(-0.5*dlgrho(1)*dz(1) + 1);
        b4  = -(i*BigOmega(1)*dz(1) + i/BigOmega(1)*dlgrho(1)*dz(1)*g);
        
        ui3 =  i*alpha_local(2)/BigOmega(2);
        ui4 = -k/BigOmega(2);

        % step 2 }}}

        name1      = 'runb.log';

        save_data1 = [mantexpnt(b1), mantexpnt(real(b2)), mantexpnt(imag(b2)), mantexpnt(b3), mantexpnt(real(b4)), mantexpnt(imag(b4))];

        name2      = 'runc.log';
        save_data2 = [mantexpnt(c1), mantexpnt(real(c2)), mantexpnt(imag(c2)), mantexpnt(c3), mantexpnt(real(c4)), mantexpnt(imag(c4))];

        name3      = 'runcw.log';
        save_data3 = [mantexpnt(cw1) mantexpnt(real(cw2)), mantexpnt(imag(cw2)), mantexpnt(real(cw3)), mantexpnt(imag(cw3)), mantexpnt(real(cw4)), mantexpnt(imag(cw4)), mantexpnt(real(cw5)), mantexpnt(imag(cw5))];

        name4      = 'runcu.log';
        save_data4 = [mantexpnt(real(cu1)), mantexpnt(imag(cu1)), mantexpnt(real(cu2)), mantexpnt(imag(cu2)), mantexpnt(real(cu3)), mantexpnt(imag(cu3)), mantexpnt(real(cu4)), mantexpnt(imag(cu4)), mantexpnt(cu5')];

        name5      = 'runcp.log';
        save_data5 = [mantexpnt(cp1), mantexpnt(real(cp2)), mantexpnt(imag(cp2)), mantexpnt(cp3), mantexpnt(cp4)];

        save(name1, 'save_data1', '-ascii');
        save(name2, 'save_data2', '-ascii');
        save(name3, 'save_data3', '-ascii');
        save(name4, 'save_data4', '-ascii');
        save(name5, 'save_data5', '-ascii');

        % initial conditions }}}

        % on construit le vecteur B des seconds membres (SM), avec :
        % - une partie conditions initiales
        % - TODO eventuellement, une partie de propagation non visqueuse (SM = 0)
        % - une partie de propagation visqueuse avec les SM appropries

        % {{{ remains

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

        sm = [sm_u'; sm_w'; sm_p'];
        sm = sm(:);

        name6      = 'runsm.log';
        save_data6 = [real(mantexpnt(sm_u)), imag(mantexpnt(sm_u)), mantexpnt(real(sm_w)), mantexpnt(imag(sm_w)), real(mantexpnt(sm_p)), imag(mantexpnt(sm_p))];
        save(name6, 'save_data6', '-ascii');
        disp 'coeff saved'
        pause

        % remains }}}

        B = [0; w(1); p(1); ... % conditions initiales
             0;    0;    0; ... % le premier pas est non visqueux
                 sm         ... % la propagation visqueuse implique des termes non nuls,
            ];                  % tous calcules a priori avec les donnees du probleme

        % construction du second membre }}}

        % {{{ viscuous K

        % helpers

        dimvisc = size(cu1, 1);
        myones  = ones( 1, dimvisc);
        myzeros = zeros(1, dimvisc);

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

        % viscuous K }}}
        
        % {{{ initial conditions
        
        ci0 = [1 ui1 ui2; ...
               0   1   0; ...
               0   0   1];

        ci1 = [0   0   0  1  ui3  ui4; ...
               0  b1  b2  0    1    0; ...
               0  b4  b3  0    0    1];

        ci  = [ci0 zeros(3, n_mod3-3) ; ci1 zeros(3, n_mod3-6)];

        % initial conditions }}}

        % {{{ sparse matrix

        preK = spdiags([diag1 diag2 diag3 diag4 diag5 diag6 diag7 diag8], ...
                       -1:6, n_mod3, n_mod3);

        size(diag1,1)
       
        preK    = preK(1:n_mod3 - 6, :);
        K       = [ci; preK];
        
        % sparse matrix }}}
 
        % {{{ solving it

        U = K\B;

        % solving it }}}

        % {{{ final processing

        % get results rows 
        u = U(1:3:end); 
        w = U(2:3:end);
        p = U(3:3:end);
        
        kz2 = (k*k).*((-g/omega/omega).*dlgrho-1) - 0.25/6.4e7;
        rho1 = (-i.*dlgrho.*w)./BigOmega;

        % final processing }}}

        % {{{ post-processing for current step

        % do some...
        clear w u p rho1 kz2;

        % post-processing for current step }}}

    end
end     % }}}

% {{{ post-processing

% do some...

% post-processing }}}

