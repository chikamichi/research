% this function processes the data sets and outputs series of plots
% which can be gathered in a movie, for example

function f = plotit()
  % advection scheme type
  scheme = 'b';

  % just to get min/max and the like
  path = ['data/' scheme '/'];
  load([path 'd_0001']);
  d_initial = d_0001;
  [imax, kmax] = size(d_initial);
  d_min = min(min(d_initial));
  d_max = max(max(d_initial));

  iso = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];

  % actually processing files
  d = dir([path 'd_*']);
  n = length(d);

  disp('computing data sets as plots...')
  progressbar
  for k = 1:n
    %disp(k)
    fname = d(k).name;
    data  = load([path fname]);
    %disp(size(data'))

    h1 = figure('visible','off');
    surf(1:imax, 1:kmax, data');
    zlim([d_min-0.05 d_max-0.5]);

    h2 = figure('visible', 'off');
    [C, hc] = contourf(1:imax, 1:kmax, data', iso);
    %text_handle = clabel(C, hc, ...);
    colorbar;

    saveas(h1, ['output/3D/'  fname '.png']);
    saveas(h2, ['output/map/' fname '.png']);

    close(h1);
    close(h2);
    clear h1;
    clear h2;
    %disp('')

    progressbar(k/n)
  end

  % get back to the shell
  quit force
end
