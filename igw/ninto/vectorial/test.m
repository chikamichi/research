clear all;
close all;

load modele1
load wind_M


res  = 1        % 10, 100, etc.
binf = 1
bsup = 9981
Z       = modele1(binf:res:bsup,1);
rho     = modele1(binf:res:bsup,2);
drho    = modele1(binf:res:bsup,3);
dlgrho  = modele1(binf:res:bsup,4);

xx         = size(Z);       % 9999
n_mod      = xx(1)          % 9999 le nombre de modes (autant que de tranches
                            % d'altitudes ?)
I(1:n_mod) = 1;

% calcul de la dérivée spatiale verticale
% TODO à vectoriser ^^
for j = 1:n_mod-1
  dz(j) = (Z(j+1)-Z(j))*1.e3;     % pas spatial en mètres, utilisé dans les
end                               % dérivées en dz j'imagine
dz(n_mod) = dz(n_mod-1);          % FIXME ajouté pour la version vectorisée, bof

dz = dz';

% vents méridiens
u0 = wind_M(:,2);
% interpolation linéaire au pas spatial de modele1
res = [];
size(u0)
u0(1:3)
u0(end-3:end)
%for i = 2:size(u0, 1)
    res = linspace(u0(1), u0(end), 21*499-499);
    %res = [res tmp];
    %tmp = [];
%end
% FIXME 500 points de trop !
size(res)
u0 = res;
u0(1:3)
u0(end-3:end)
%size(u0)
%size(dz)
%u0diff = diff(u0);
%size(u0diff)
%alpha = diff(u0)./dz;
%size(alpha)

