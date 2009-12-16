% this function processes the data sets and outputs series of plots
% which can be gathered in a movie, for example

function f = plotit()
  % 3D plots mesh resolution
  dmesh    = 1;
  % available data sampling resolution
  % (should be the same as mod_backup in td4.f90)
  timestep = 10;

  % advection scheme type
  scheme = 'b';
  path = ['data/' scheme '/'];

  % get the files list
  d = dir([path 'd_*']);
  n = length(d);

  % bold assumption: get some min/max ranges using the last computation
  %d_final = d(n).name;           % beware, it's 2D actually
  % FIXME bold fix cuz I dunno how to get the var name expanded within the statements
  [imax, kmax] = size(load([path d(n).name]));
  zmax = 2*n*timestep;
  %zmax  = 2;
  %d_min = min(min(load([path d(n).name])));
  %d_max = max(max(load([path d(n).name])));
  %means = load('mean.dat');

  iso = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 2, 5, 7, 10, 20, 30, 50, 100, 200, 500];
  %iso = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5];
  % build a logarithmic tracer concentrations scale on purpose
  %scaling_low  = 0.8;     % this is to give some latitude to the log scale
  %scaling_high = 1.2;     % so the borders are not the min and max log values
  %n_log_values = 30;      % let's get n_log_values out of it
  %if d_min == 0
    %d_min = 1e-3;
  %end
  %iso = linspace(log(d_min*scaling_low), log(d_max*scaling_high), n_log_values);

  % actually processing files
  disp('computing data sets as plots...')
  progressbar
  for k = 1:n
    %disp(k)
    fname = d(k).name;
    data  = load([path fname]);
    %disp(size(data'))

    h1 = figure('visible','off');
    surf(1:dmesh:imax, 1:dmesh:kmax, data(1:dmesh:end, 1:dmesh:end)');
    axis([0 imax 0 kmax 0 zmax]); % constant axis range
    %caxis([0, zmax]); % btw colorbar range must be constant for the 
                        % shading interpolation to remain consistent
                        % along the movie
                        % ARGH: actually wrong! that's not iso here!
    zlim([0 zmax]);
    shading interp;
    grid on;
    view([110 40]);
    colorbar;

    %h2 = figure('visible', 'off');
    %%[C, hc] = contourf(1:imax, 1:kmax, data');
    %[C, hc] = contourf(1:imax, 1:kmax, data', iso);
    %% if one want to put labels on the iso:
    %%text_handle = clabel(C, hc, ...);
    %%set(gca, 'zscale', 'log');
    %colorbar;

    saveas(h1, ['output/3D/'  fname '.png']);
    %saveas(h2, ['output/map/' fname, '.png']);

    %%mean_fig = contourf(1:imax, 1:kmax, means', iso);
    %%colorbar;
    %%saveas(mean_fig, ['means_contours.png']);

    close(h1);
    %close(h2);
    %%close(mean_fig);
    clear h1;
    %clear h2;
    %%clear mean_fig;
    %%disp('')

    progressbar(k/n)
  end

  % get back to the shell
  quit force
end
