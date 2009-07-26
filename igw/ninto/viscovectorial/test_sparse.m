% script de test pour la construction d'une matrice creuse

% dimensions

n_mod  = 5;
n_mod3 = 3*n_mod;

% coefficients et pré-diagonales

u1   = 1;
p1   = 1;
ui1  = 3;
ui2  = 4;
ui3  = 5;
ui4  = 6;
b1   = 7;
b2   = 8;
b3   = 9;
b4   = 10;
%c1   = repmat(11,     n_mod,1);
%c2   = repmat(12,     n_mod,1);
%c3   = repmat(13,     n_mod,1);
%c4   = repmat(14,     n_mod,1);
cu1  = repmat(15,    n_mod,1);
cu2  = repmat(16,    n_mod,1);
cu3  = repmat(17,    n_mod,1);
cu4  = repmat(18,    n_mod,1);
cu5  = repmat(19,    n_mod,1);
cw1  = repmat(20,    n_mod,1);
cw2  = repmat(21,    n_mod,1);
cw3  = repmat(22,    n_mod,1);
cw4  = repmat(23,    n_mod,1);
cw5  = repmat(24,    n_mod,1);
cp1  = repmat(25,    n_mod,1);
cp2  = repmat(26,    n_mod,1);
cp3  = repmat(27,    n_mod,1);
cp4  = repmat(28,    n_mod,1);

% helpers

dimvisc = size(cu1, 1);
myones  = ones(1, dimvisc);
myzeros = zeros(1,dimvisc);

% diagonales alternées (les * dénotent celles affectées d'une permutation circulaire)

diag1   = [cw5';  cp4'; myzeros];       % *
diag1   = diag1(:);

diag2   = [cu2'; -myones; -myones];
diag2   = diag2(:);

diag3   = [myzeros myzeros myzeros];    % *
diag3   = diag3';

diag4   = [cw4'; cp2'; myzeros];        %*
diag4   = diag4(:);

diag5   = [cu1'; cw1'; cp1'];
diag5   = diag5(:);

diag6   = [myzeros; cu5'; cw2'];        % *
diag6   = diag6(:);

diag7   = [cw3'; cp3'; cu4'];         % *
diag7   = diag7(:);

diag8   = [cu3'; myones; myones];
diag8   = diag8(:);

% matrices des conditions initiales

ci0 = [1 ui1 ui2; ...
       0   1   0; ...
       0   0   1];

ci1 = [0   0   0  1  ui3  ui4; ...
       0  b1  b2  0    1    0; ...
       0  b4  b3  0    0    1];

ci  = [ci0 zeros(3, n_mod3-3) ; ci1 zeros(3, n_mod3-6)];

% construction de K

preK = spdiags([diag1 diag2 diag3 diag4 diag5 diag6 diag7 diag8], ...
               -1:6, n_mod3, n_mod3);
%full(preK)

% troncature de la base, remplacée par les conditions initiales (6 lignes)

preK = preK(1:n_mod3 - 6, :);
K    = [ci; preK];
%disp('K :')
%size(K)

%figure
%spy(K)
full(K)

