% TEST GRAVITY WAVES %

clear all;
close all;

% %load /Users/ninto/occhipinti/MODELE/USSA76/rho_drho_dlgrho
% %load /Users/ninto/occhipinti/MODELE/USSA76/RHOout
% load /Users/ninto/occhipinti/MODELE/USSA76/modele1
load Venus_atm_jour.dat

% EARTH
% %Z=RHOout(:,1);
% Z=modele1(:,1);
% VENUS
Z= interp(Venus_atm_jour(31:71,1),80);

nnn=max(size(Z));

% EARTH
% %rho=1.225*exp(-Z*1.e3/8.e3);
% %rho=RHOout(:,2);
% rho=modele1(:,2);
% VENUS
rho=interp(Venus_atm_jour(31:71,3),80);

% EARTH
% %drho=-1/8.e3*rho;   %gradient(rho,(Z(2)-Z(1))*1.e3);      
% %drho=RHOout(:,3);
% drho=modele1(:,3);
% VENUS
drho=gradient(rho,(Z(2)-Z(1))*1.e3);

% EARTH
% %dlgrho(1:nnn)=-1/8.e3; %gradient(log(rho),(Z(2)-Z(1))*1.e3);    
% %dlgrho=RHOout(:,4);
% dlgrho=modele1(:,4);
% VENUS
dlgrho=gradient(log(rho),(Z(2)-Z(1))*1.e3); 

% VENUS WINDS
% % zonal wind
%wind=-1*spline(Venus_atm_jour(31:71,1),Venus_atm_jour(31:71,13),Z);
% alpha = 0;
% merional wind
wind=-1*spline(Venus_atm_jour(31:71,1),Venus_atm_jour(31:71,14),Z);
grad_wind = gradient(wind,245) ;
%alpha = 0; % between 31 and 48
%alpha = -2.89e-3; % between 49 and 71

xx=size(Z);
n_mod=xx(1);
I(1:n_mod)=1;

for j=1:n_mod-1
    dz(j)=(Z(j+1)-Z(j))*1.e3;
end

% CONSTANT GRADIANT WIND
% u0b=0;
% u0t=1;
% alpha=(u0t-u0b)/(max(Z)-min(Z)*1.e3);

M=5.9768e+24;
G=6.67e-11;
R=6378.e3;

%----------------------------------
%   omega, kx & ky from file
%----------------------------------
% 
% load /Users/ninto/occhipinti/SRC/TGW/RUN/grid
% 
%            %TIME SPACE%
% lon=[grid(1,1):grid(1,3):grid(1,2)];
% lat=[grid(2,1):grid(2,3):grid(2,2)];
% lat2=sort(lat,'descend');
% time=[grid(3,1):grid(3,3):grid(3,2)];
% 
% clear grid;
% 
%          %FREQUENCY SPACE%
% T=max(time)-min(time);
% dt=time(2)-time(1);
% wn=[-1/2/dt:1/T:1/2/dt]*2*pi;
% 
% deg2m=40.e6/360;
% 
% Lx=(max(lon)-min(lon))*deg2m;
% dlx=abs(lon(2)-lon(1))*deg2m;
% kx=[-1/2/dlx:1/Lx:1/2/dlx]*2*pi;
% 
% Ly=(max(lat)-min(lat))*deg2m;
% dly=abs(lat(2)-lat(1))*deg2m;
% ky=[-1/2/dly:1/Ly:1/2/dly]*2*pi;
% %----------------------------------

%----------------------------------
%   omega, kx & ky à la main
%----------------------------------

% EARTH & TSUNAMI
% %T=[1:10:30]*60
% T=[50 100 150]*60
% wn=2*pi./T;
% 
% h=2500;
% g=9.8;
% c=sqrt(g*h);
% kx=wn/c;
% ky=kx;

% VENUS
T=[5:5:30]*60
%T=[10]*60; % 0 100 150
wn=2*pi./T;

L = [50 100 150 200 250 300 350 400]*1.e3;
kx = 2*pi./L;
ky=kx;
h=2500;
% g=9.8;
% c=sqrt(g*h);
% kx=wn/c;
% ky=kx;


%----------------------------------

nlat=max(size(ky));
nlon=max(size(kx));
nt=max(size(wn));

%----------------------------------
% 
%----------------------------------

load BestView2
%----------------------------------
for ik=1:nlon %(nlon-1)/2+2:nlon
    ik
    drawnow;

    k=kx(ik);

    %      % FIGURE 3D
    %      figure(2*ik);
    %      set(2*ik,'position',[44   250   500   750]);
    %      set(gcf,'Color','w');
    %      set(gca,'LineWidth',[2])
    %      set(gca,'fontsize',18);
    %      view(VV);
    %      title(['Hwater=',num2str(h),'m and kx=', num2str(k),' rad/m'],'fontsize',18);
    % %     xlabel('omega w (rad/s)','fontsize',18);
    %      xlabel('period T (min)','fontsize',18);     
    %      ylabel('Vertical V_r_e_a_l (m/s)','fontsize',18);
    %      Zlabel('altitude (km)','fontsize',18);
    %      grid on;
    %      hold on;
    % 
    %      figure(2*ik+1);
    %      set(2*ik+1,'position',[44   250   500   750]); 
    %      set(gcf,'Color','w');
    %      set(gca,'LineWidth',[2])
    %      set(gca,'fontsize',18);
    %      view(VV);    
    %      title(['Hwater=',num2str(h),'m   and   kx=', num2str(k),' rad/m'],'fontsize',18);
    % %    xlabel('omega (rad/s)','fontsize',18);
    %      xlabel('period T (min)','fontsize',18);
    %      ylabel('Vertical V_i_m_a_g (m/s)','fontsize',18);
    %      Zlabel('altitude (km)','fontsize',18);
    %      grid on;
    %      hold on;     


    for iw=1:nt %(nt-1)/2+2:nt

        %%ik=(nlon-1)/2+1 + 20
        %%iw=(nt-1)/2+1 + 20

        w(1)=1;
        % p(1)=-1;

        %ik-(nlon-1)/2-1
        % k=kx(ik);
        % k=2*pi/(23.8095*60*sqrt(9.8*2.e3));

        %iw-(nt-1)/2-1
        omega=wn(iw);
        % omega=2*pi/(23.8095*60); 

        g=9.8; %G*M/(R+Z(10)*1.e3)^2;

        cz(ik,iw,j)=omega/k/sqrt(-g/omega/omega*dlgrho(j)-1);


        for j=1:1 % n_mod-1

            %     % CONSTANT GRADIENT WIND
            %     u0(j)=Z(j)*alpha+u0b;

            u0(j)= wind(j);
            alpha = grad_wind(j);
            %     if Z >= 94.001 
            %         alpha = 3.e-4;
            %     end

            g=9.8; %G*M/(R+Z(j)*1.e3)^2;


            kz2(j)=k*k*(-g/omega/omega*dlgrho(j)-1) -0.25/6.4e7;
            N2(j)=-g*dlgrho(j);

            %   eq. 4*
            rho1(j)=-i/(omega-k*u0(j))*dlgrho(j)*w(j);

            if N2(j) >= omega*omega 

                %   cas propagative 
                %   eq. P
                p(j)=i/(k*k)*(alpha*k-i*sqrt(kz2(j))*(omega-k*u0(j))-0.5*(omega-k*u0(j))*dlgrho(j))*w(j);

            else

                %   cas non propagative 
                %   eq. P
                p(j)=i/(k*k)*(alpha*k-sqrt(kz2(j))*(omega-k*u0(j))-0.5*(omega-k*u0(j))*dlgrho(j))*w(j);

            end

            %   eq. 1*
            u(j)=1/(omega-k*u0(j))*(-i*w(j)*alpha+k*p(j));


            %   eq. 3*
            w(j+1)=-i*k*u(j)*dz(j)+0.5*dlgrho(j)*dz(j)*w(j)+w(j);
            %    w(j+1)=w(j);

            %   eq. 2*    
            p(j+1)=i*(omega-u0(j)*k)*dz(j)*w(j)-0.5*dlgrho(j)*dz(j)*p(j)-rho1(j)*dz(j)*g+p(j);    
            %    p(j+1)=p(j);

        end


        for j=2:n_mod-1

            %     % CONSTANT GRADIENT WIND
            %     u0(j)=Z(j)*alpha+u0b;

            u0(j)= wind(j);
            alpha = grad_wind(j);
            %     if Z >= 94.001 
            %         alpha = 3.e-4;
            %     end

            g=9.8; % G*M/(R+Z(j-1)*1.e3)^2;

            %k2(j)=k*k;
            kz2(j)=k*k*(-g/omega/omega*dlgrho(j)-1) -0.25/6.4e7;
            N2(j)=-g*dlgrho(j);

            %   eq. 1*
            u(j)=1/(omega-k*u0(j))*(-i*alpha*w(j)+k*p(j));

            %   eq. 4*
            rho1(j)=-i/(omega-k*u0(j))*dlgrho(j)*w(j);

            %   eq. 3*
            w(j+1)=-i*k*u(j)*2*dz(j)+0.5*dlgrho(j)*2*dz(j)*w(j)+w(j-1);   %+4*w(j);

            %   eq. 2*    
            p(j+1)=i*(omega-u0(j)*k)*2*dz(j)*w(j)-0.5*dlgrho(j)*2*dz(j)*p(j)-rho1(j)*2*dz(j)*g+p(j-1); %+4*p(j);    

        end


        kz2(n_mod)=k*k*(-g/omega/omega*dlgrho(n_mod)-1) -0.25/6.4e7; 
        N2(n_mod)=-g*dlgrho(n_mod);

        u0(n_mod)= wind(n_mod);
        %u0(n_mod)=Z(n_mod)*alpha+u0b;

        % eq. 1*
        u(n_mod)=1/(omega-k*u0(n_mod))*(-i*alpha*w(n_mod)+k*p(n_mod));

        % eq. 4*
        rho1(n_mod)=-i/(omega-k*u0(n_mod))*dlgrho(n_mod)*w(n_mod);

        %------------------------------------------------------------
        %               da usare con file grid
        %
        % wmax(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=max(abs(real(w)));
        % umax(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=max(abs(real(u)));
        % pmax(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=max(abs(real(p)));
        % rhomax(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=max(abs(real(rho1)));
        % KZ(ik-(nlon-1)/2-1,iw-(nt-1)/2-1,:)=sqrt(kz2);
        % 
        % if max(abs(real(w))) >= 2.
        %     ctrw(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=NaN;
        % else
        %     ctrw(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=1.;
        % end
        % 
        % if sqrt(kz2(j)) >= 2*pi/50.e3 
        % 
        %     A(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=1;
        %     
        % else
        % 
        %     A(ik-(nlon-1)/2-1,iw-(nt-1)/2-1)=NaN;
        %     
        % end
        %------------------------------------------------------------

        %------------------------------------------------------------
        %            da usare con kx e omega ˆ la main
        %
        wmax(ik,iw)=max(abs(real(w)));
        umax(ik,iw)=max(abs(real(u)));
        pmax(ik,iw)=max(abs(real(p)));
        rhomax(ik,iw)=max(abs(real(rho1)));
        KZ(ik,iw,:)=sqrt(kz2);

        if max(abs(real(w))) >= 2.
            ctrw(ik,iw)=NaN;
        else
            ctrw(ik,iw)=1.;
        end

        if sqrt(kz2(j)) >= 2*pi/50.e3 

            A(ik,iw)=1;

        else

            A(ik,iw)=NaN;

        end

        %------------------------------------------------------------

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure(100);
        set(100,'position',[10 1 750 600])
        set(gcf,'Color','w');
        fnt = 20;
        fnt2 = 18;
        Zmax = 140;
        Zmin = 60;
        subplot(1,4,1);
        plot(real(w),Z,'r','linewidth',2);
        hold on;
        plot(imag(w),Z,'g','linewidth',2);
        hold off;
        v=axis;
        axis([v(1) v(2) Zmin Zmax])
        %axis([-2 2 Zmin Zmax])
        xlabel('Vertical V (m/s)','fontsize',fnt);
        ylabel('Altitude (km)','fontsize',fnt);
        %title(['omega= ',num2str(omega),' rad/s'],'fontsize',fnt);
        title(['T= ',num2str(2*pi/omega/60),' min'],'fontsize',fnt);
        set(gca,'fontsize',fnt2);
        subplot(1,4,2);
        plot(real(u),Z,'r','linewidth',2);
        hold on;
        plot(imag(u),Z,'g','linewidth',2);
        hold off;
        v=axis;
        axis([v(1) v(2) Zmin Zmax])
        xlabel('Horizontal V (m/s)','fontsize',fnt);
        set(gca,'fontsize',fnt2);
        subplot(1,4,3);
        plot(real(p*1.e-3),Z,'r','linewidth',2);
        hold on;
        plot(imag(p*1.e-3),Z,'g','linewidth',2);
        hold off;
        v=axis;
        axis([v(1) v(2) Zmin Zmax])
        %axis([-6000 6000 Zmin Zmax])
        xlabel('Pressure (kPa)','fontsize',fnt);
        set(gca,'fontsize',fnt2);
        subplot(1,4,4);
        plot(real(rho1),Z,'r','linewidth',2);
        hold on;
        plot(imag(rho1),Z,'g','linewidth',2);
        hold off;
        v=axis;
        axis([v(1) v(2) Zmin Zmax])
        %axis([-0.2 0.2 Zmin Zmax])
        xlabel('Density (kg/m3)','fontsize',fnt);
        %title(['kx=', num2str(k),' rad/m'],'fontsize',fnt);
        title(['\lambda_x=', num2str(fix(2*pi/k*1.e-3*10)/10),' km'],'fontsize',fnt);
        set(gca,'fontsize',fnt2);

        figure(100);
        F=getframe(gcf);
        [X,Map] = frame2im(F);
        file = ['ModeVenus_L' num2str(round(1/k*1.e-3)) '_T_' num2str(2*pi/omega/60) '.tif'];
        imwrite(X,file,'tif','Compression','none');

        % T(iw)/60
        % L(ik)*1.e-3
        % disp('PRESS ENTER TO CONTINUE...');
        %pause
        close

        Aw_max_real(ik,iw) = max(abs(real(w(1:1201))));
        Aw_max_imag(ik,iw) = max(abs(imag(w(1:1201))));

        Aw_max_real2(ik,iw) = max(abs(real(w(2401:3201))));
        Aw_max_imag2(ik,iw) = max(abs(imag(w(2401:3201))));

        Aw_max_real0(ik,iw) = abs(real(w(1)));
        Aw_max_imag0(ik,iw) = abs(imag(w(1)));


        Au_max_real(ik,iw) = max(abs(real(u(1:1201))));
        Au_max_imag(ik,iw) = max(abs(imag(u(1:1201))));

        Au_max_real2(ik,iw) = max(abs(real(u(2401:3201))));
        Au_max_imag2(ik,iw) = max(abs(imag(u(2401:3201))));

        Au_max_real0(ik,iw) = abs(real(u(1)));
        Au_max_imag0(ik,iw) = abs(imag(u(1)));


        Ap_max_real(ik,iw) = max(abs(real(p(1:1201)*1.e-3)));
        Ap_max_imag(ik,iw) = max(abs(imag(p(1:1201)*1.e-3)));

        Ap_max_real2(ik,iw) = max(abs(real(p(2401:3201)*1.e-3)));
        Ap_max_imag2(ik,iw) = max(abs(imag(p(2401:3201)*1.e-3)));

        Ap_max_real0(ik,iw) = abs(real(p(1)*1.e-3));
        Ap_max_imag0(ik,iw) = abs(imag(p(1)*1.e-3));


        Arho_max_real(ik,iw) = max(abs(real(rho1(1:1201))));
        Arho_max_imag(ik,iw) = max(abs(imag(rho1(1:1201))));

        Arho_max_real2(ik,iw) = max(abs(real(rho1(2401:3201))));
        Arho_max_imag2(ik,iw) = max(abs(imag(rho1(2401:3201))));

        Arho_max_real0(ik,iw) = abs(real(rho1(1)));
        Arho_max_imag0(ik,iw) = abs(imag(rho1(1)));

        if (Aw_max_real2(ik,iw) >= 500 | Aw_max_imag2(ik,iw) >= 500)
            disp('top --- >>> NaN')

            Aw_max_real2(ik,iw)     = NaN;    
            Au_max_real2(ik,iw)     = NaN;
            Ap_max_real2(ik,iw)     = NaN;
            Arho_max_real2(ik,iw)   = NaN;

        end

        if (Aw_max_real(ik,iw) >= 50 | Aw_max_imag(ik,iw) >= 50)
            disp('low --- >>> NaN')

            Aw_max_real(ik,iw)     = NaN;    
            Au_max_real(ik,iw)     = NaN;
            Ap_max_real(ik,iw)     = NaN;
            Arho_max_real(ik,iw)   = NaN;

        end


        % FIGURE INFO
        figure(300);
        subplot(1,3,1);
        plot(real(sqrt(kz2)),Z,'r');
        hold on;
        plot(imag(sqrt(kz2)),Z,'g');
        xlabel('kz (rad/m)','fontsize',12);
        ylabel('Altitude (km)','fontsize',16);
        hold off;
        subplot(1,3,2);
        plot(sqrt(N2),Z,'r');
        hold on;
        v=axis;
        plot([omega omega],[v(3) v(4)],'g');
        hold off;
        xlabel('N and omega (rad/s)','fontsize',12);
        subplot(1,3,3);
        plot(2*pi./real(sqrt(kz2))*1.e-3,Z,'r');
        hold on;
        plot(2*pi./imag(sqrt(kz2))*1.e-3,Z,'g');
        hold off;
        xlabel('lambda Z (km)','fontsize',12);
        
        disp('PRESS ENTER TO CONTINUE...');
        pause
        close


        figure(200);
        subplot(1,3,1);
        plot(real(w),Z,'r');
        hold on;
        plot(imag(w),Z,'g');
        hold off;
        v2=axis;
        xlabel('Vertical V Num (m/s)','fontsize',12);
        ylabel('Altitude (km)','fontsize',16);
        title(['omega= ',num2str(omega),' rad/s'],'fontsize',12);
        subplot(1,3,2);
        plot(imag(exp(-i*sqrt(kz2(1))*Z*1.e3)),Z,'g');
        hold on;
        plot(real(exp(-i*sqrt(kz2(1))*Z*1.e3)),Z,'r');
        hold off;
        axis(v2);
        xlabel('Vertical V theo (m/s)','fontsize',12);
        subplot(1,3,3);
        plot(real(w)-real(exp(-i*sqrt(kz2(1))*Z*1.e3))',Z,'r');
        hold on;
        plot(imag(w)-imag(exp(-i*sqrt(kz2(1))*Z*1.e3))',Z,'g');
        hold off;
        %axis(v2);
        xlabel('Vertical V res (m/s)','fontsize',12);
        title(['kx=', num2str(k),' rad/m'],'fontsize',12);
        
        pause
        close
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % FIGURE 3D
        figure(2*ik);
        plot3(I*2*pi/omega/60,real(w),Z,'k');
        plot3(I*2*pi/omega/60,imag(w),Z,'r');
        title(['Hwater=',num2str(h),'m and kx=', num2str(k),' rad/m'],'fontsize',18);
        xlabel('omega (rad/s)','fontsize',18);
        ylabel('Vertical V (m/s)','fontsize',18);
        Zlabel('altitude (km)','fontsize',18);
        % hold on;
        
        figure(2*ik+1);
        plot3(I*2*pi/omega/60,imag(w),Z,'k');
        title(['Hwater=',num2str(h),'m   and   kx=', num2str(k),' rad/m'],'fontsize',18);
        xlabel('omega (rad/s)','fontsize',18);
        ylabel('Vertical V (m/s)','fontsize',18);
        Zlabel('altitude (km)','fontsize',18);
        hold on;
        % END FIGURE 3D

        xxx=input('...clicca !!!');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clear w u p rho1 kz2

    end 
end

Nsol=sqrt(N2(1));
Nmax=max(sqrt(N2));
Nmin=min(sqrt(N2));

njump=0; 
njump2=size(wmax);

% n1=(nt-1)/2+2+njump
% n2=(nlon-1)/2+2

n1=1;
n2=1;

h=figure;
set(h,'position',[2 700 800 600]);
set(gcf,'Color','w');
font_size = 22;
font_sizeN = 18;

flag_plot = input('-Plot omega-k or T-lambda ? (0 or 1) '); 

subplot(2,2,1);
if flag_plot == 0
    pcolor(wn(n1:nt),kx(n2:nlon),log10(wmax(:,1+njump:njump2(2))));
else
    pcolor(2*pi./wn(n1:nt)/60,2*pi./kx(n2:nlon)*1.e-3,log10(wmax(:,1+njump:njump2(2))));
end
shading flat;
v=axis;
hold on;
plot([Nmax Nmax],[v(3) v(4)],'k');
plot([Nsol Nsol],[v(3) v(4)],'k');
plot([Nmin Nmin],[v(3) v(4)],'k');
colorbar;
title('W_r_e_a_l','fontsize',font_size);
if flag_plot == 0
    ylabel('kx (rad/m)','fontsize',font_size);
else
    ylabel('lambda (km)','fontsize',font_size);
end
set(gca,'LineWidth',[2])
set(gca,'fontsize',font_sizeN);

subplot(2,2,2);
if flag_plot == 0
    pcolor(wn(n1:nt),kx(n2:nlon),log10(umax(:,1+njump:njump2(2))));
else
    pcolor(2*pi./wn(n1:nt)/60,2*pi./kx(n2:nlon)*1.e-3,log10(umax(:,1+njump:njump2(2))));
end
shading flat;
v=axis;
hold on;
plot([Nmax Nmax],[v(3) v(4)],'k');
plot([Nsol Nsol],[v(3) v(4)],'k');
plot([Nmin Nmin],[v(3) v(4)],'k');
colorbar;
title('U_r_e_a_l','fontsize',font_size);
set(gca,'LineWidth',[2])
set(gca,'fontsize',font_sizeN);

subplot(2,2,3);
if flag_plot == 0
    pcolor(wn(n1:nt),kx(n2:nlon),log10(pmax(:,1+njump:njump2(2))));
else
    pcolor(2*pi./wn(n1:nt)/60,2*pi./kx(n2:nlon)*1.e-3,log10(pmax(:,1+njump:njump2(2))));
end
shading flat;
v=axis;
hold on;
plot([Nmax Nmax],[v(3) v(4)],'k');
plot([Nsol Nsol],[v(3) v(4)],'k');
plot([Nmin Nmin],[v(3) v(4)],'k');
colorbar;
if flag_plot == 0
    ylabel('kx (rad/m)','fontsize',font_size);
    xlabel('omega (rad/s)','fontsize',font_size);
else
    ylabel('lambda (km)','fontsize',font_size);
    xlabel('period T (min)','fontsize',font_size);
end
title('P_r_e_a_l','fontsize',font_size);
set(gca,'LineWidth',[2])
set(gca,'fontsize',font_sizeN);

subplot(2,2,4);
if flag_plot == 0
    pcolor(wn(n1:nt),kx(n2:nlon),log10(rhomax(:,1+njump:njump2(2))));
else
    pcolor(2*pi./wn(n1:nt)/60,2*pi./kx(n2:nlon)*1.e-3,log10(rhomax(:,1+njump:njump2(2))));
end
shading flat;
v=axis;
hold on;
plot([Nmax Nmax],[v(3) v(4)],'k');
plot([Nsol Nsol],[v(3) v(4)],'k');
plot([Nmin Nmin],[v(3) v(4)],'k');
colorbar;
if flag_plot == 0
    xlabel('omega (rad/s)','fontsize',font_size);
else
    xlabel('period T (min)','fontsize',font_size);
end
title('Rho_r_e_a_l','fontsize',font_size);
set(gca,'LineWidth',[2])
set(gca,'fontsize',font_sizeN);



% FIGURE EMELINE

% %Aw_max_real(2,1)        = NaN;
% Aw_max_real(3,1)        = NaN;
% Aw_max_real(4,[1])      = NaN;  % [1 2]
% Aw_max_real(5,[1 2])    = NaN;  % []
% Aw_max_real(6,[1 2])    = NaN;  % [1 2 3]
% Aw_max_real(7,[1 2])    = NaN;  % [1 2 3]
% Aw_max_real(8,[1 2 3])  = NaN;  % [1 2 3 4]
% 
% %Aw_max_real2(2,1)        = NaN;
% Aw_max_real2(3,1)       = NaN;
% Aw_max_real2(4,[1])     = NaN;
% Aw_max_real2(5,[1 2])   = NaN;
% Aw_max_real2(6,[1 2])   = NaN;
% Aw_max_real2(7,[1 2])   = NaN;
% Aw_max_real2(8,[1 2 3]) = NaN;
% 
% %Au_max_real(2,1)        = NaN;
% Au_max_real(3,1)        = NaN;
% Au_max_real(4,[1])      = NaN;
% Au_max_real(5,[1 2])    = NaN;
% Au_max_real(6,[1 2])    = NaN;
% Au_max_real(7,[1 2])    = NaN;
% Au_max_real(8,[1 2 3])  = NaN;
% 
% %Au_max_real2(2,1)        = NaN;
% Au_max_real2(3,1)       = NaN;
% Au_max_real2(4,[1])     = NaN;
% Au_max_real2(5,[1 2])   = NaN;
% Au_max_real2(6,[1 2])   = NaN;
% Au_max_real2(7,[1 2])   = NaN;
% Au_max_real2(8,[1 2 3]) = NaN;

load COLOR.mat

% Vitesse Verticale
figure2
colormap(COLOR)
font_size = 24;
subplot(2,2,1)
%imagesc2(T/60,L*1.e-3,log10(Aw_max_real./Aw_max_real0));
%%pcolor(T/60,L*1.e-3,log10(Aw_max_real2./Aw_max_real0));
pcolor(T/60,L*1.e-3,Aw_max_real2./Aw_max_real0);
caxis([1 3])
%caxis([0 2.5])
colorbar
set(gca,'fontsize',18)
set(gca,'LineWidth',[2])
title('V_z amplification upper atmo.','fontsize',font_size);
xlabel('period T (min)','fontsize',font_size);
ylabel('lambda (km)','fontsize',font_size);
subplot(2,2,3)
%imagesc2(T/60,L*1.e-3,log10(Aw_max_real2./Aw_max_real0));
%%pcolor(T/60,L*1.e-3,log10(Aw_max_real./Aw_max_real0));
pcolor(T/60,L*1.e-3,Aw_max_real./Aw_max_real0);
caxis([1 3])
%caxis([0 1])
colorbar
set(gca,'fontsize',18)
set(gca,'LineWidth',[2])
title('V_z amplification low atmo.','fontsize',font_size);
xlabel('period T (min)','fontsize',font_size);
ylabel('lambda (km)','fontsize',font_size);

% Vitesse Horizontale
subplot(2,2,2)
%imagesc2(T/60,L*1.e-3,log10(Au_max_real./Au_max_real0));
%%pcolor(T/60,L*1.e-3,log10(Au_max_real2./Au_max_real0));
pcolor(T/60,L*1.e-3,Au_max_real2./Au_max_real0);
caxis([1 10])
%caxis([0 2.5])
colorbar
set(gca,'fontsize',18)
set(gca,'LineWidth',[2])
title('V_h amplification upper atmo.','fontsize',font_size);
xlabel('period T (min)','fontsize',font_size);
subplot(2,2,4)
%imagesc2(T/60,L*1.e-3,log10(Au_max_real2./Au_max_real0));
%%pcolor(T/60,L*1.e-3,log10(Au_max_real./Au_max_real0));
pcolor(T/60,L*1.e-3,Au_max_real./Au_max_real0);
caxis([1 3])
%caxis([0 1])
colorbar
set(gca,'fontsize',18)
set(gca,'LineWidth',[2])
title('V_h amplification low atmo.','fontsize',font_size);
xlabel('period T (min)','fontsize',font_size);

