% script de test pour la construction du second membre

n_mod  = 7;
n_mod3 = 3*n_mod;
visc   = 1.3e-5;
phi    = 2.*visc.*sqrt(rho);

% on construit le vecteur B des seconds membres (SM), avec :
% - une partie conditions initiales
% - TODO éventuellement, une partie de propagation non visqueuse (SM = 0)
%   ce qui nécessite de modifier K !
% - une partie de propagation visqueuse avec les SM appropriés

% calcul du second membre

% TODO il faut éventuellement calculer un linspace des u0 à partir du gradient linéaire, si on utilise pas de profil de vents. Pour les profils de vents, s'il n'y a pas le même nombre de pas en espace, il faut également faire des linspaces entre les points de donnée (interpolations linéaires)

% il y a 6 lignes de conditions initiales, donc deux par inconnues
for j = 3:n_mod-1
    sm_tmp = (phi(j)/(BigGamma(j)*dz(j)))*expniphase*(u0(j+1) - 2*u0(j) + u0(j-1));   
end
sm_tmp(n_mod) = (phi(n_mod)/(BigGamma(n_mod)*dz(n_mod-1)))*expniphase*(u0(n_mod) - 2*u0(n_mod-1) + u0(n_mod-2));
sm_tmp = sm_tmp(3:end);

sm_u   = i.*sm_tmp./2;
sm_w   = k.*sm_tmp;
sm_p   = zeros(1,n_mod-2);

% vecteur des valeurs alternées
sm = [sm_u'; sm_w'; sm_p'];
sm = sm(:);

%disp('sm :')
%size(sm)

B = [u(1); w(1); p(1); ... % conditions initiales
        0;    0;    0; ... % le premier pas est non visqueux
            sm         ... % la propagation visqueuse implique des termes non nuls,
    ];                     % tous calculés a priori avec les données du problème

%disp('B :')
%size(B)

