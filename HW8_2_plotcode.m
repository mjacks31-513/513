
% The following is a modified version of fdtd_original.m obtained from
% https://www.mathworks.com/matlabcentral/fileexchange/7459-fdtd1d-m
% The additions include the ability to extract the field at certain
% locations and to define different runs.
% To speed up the run, comment out the plotting in the loop.
% The basics of the algorithm when sigma = 0 is described in
% https://my.ece.utah.edu/~ece6340/LECTURES/lecture%2014/FDTD.pdf
% See also https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scott Hudson, WSU Tri-Cities
%1D electromagnetic finite-difference time-domain (FDTD) program.
%Assumes Ey and Hz field components propagating in the x direction.
%Fields, permittivity, permeability, and conductivity
%are functions of x. Try changing the value of "profile".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultTextFontName','Times');
set(0,'DefaultTextFontSize',16);
set(0,'DefaultAxesFontSize',16);

% https://www.mathworks.com/matlabcentral/answers/338733-how-to-stop-legend-from-adding-data1-data2-when-additional-data-is-plotted
% set(0,'DefaultLegendAutoUpdate','off');

%close all;
clear all;

animate = 0;

eps0 = 8.854e-12; % permittivity of free space
mu0  = pi*4e-7;   % permeability of free space
run = 1;
c = sqrt(1/eps0/mu0);

if run == 0
    profile = 0; % eps = eps_o, mu = mu_o, sigma = 0.
    source = 2;  % Gaussian pulse at left boundary
    xg = Lx;
    Niter = 400; % # of iterations to perform
end

if run == 1
    profile = 1;
    source  = 1;
    Lx  = 5;       % Domain length in meters
    Nx  = 500;     % Spatial samples in domain
    ixb = Nx/2;
    fs = 300e6;   % Source frequency in Hz
    lamda = Nx*c/fs/Lx;
    fstr = '300 MHz';
    Niter = 500;  % Number of iterations to perform
    ip = Nx/2 - 2; % Index of probe
    ylims = [-3, 3];
end

ds = Lx/Nx; % spatial step in meters
dt = ds/fs; % "magic time step"
% See https://my.ece.utah.edu/~simpson/ECE5340/Taflove%20Chpt.%202.pdf
% for definition of magic time step.

% Scale factors for E and H
ae = ones(Nx,1)*dt/(ds*eps0);
am = ones(Nx,1)*dt/(ds*mu0);
as = ones(Nx,1);

% Create grid of epsilon, mu, sigma.
[epsr,mur,sigma] = fdtd_profile(profile, Nx, ixb);

figure(1);
fdtd_profile_plot(profile, Nx, ixb);

ae = ae./epsr;
am = am./mur;
ae = ae./(1+dt*(sigma./epsr)/(2*eps0));
as = (1-dt*(sigma./epsr)/(2*eps0))./(1+dt*(sigma./epsr)/(2*eps0));

% Initialize fields to zero
Hz = zeros(Nx,1);
Ey = zeros(Nx,1);
Ey_forward = zeros(ixb-1, 1);
Hz_forward = zeros(ixb-1, 1);
figure(2);clf
    set(gcf,'doublebuffer','on'); % For smoother graphics
    grid on;
    plot(Ey,'b','LineWidth',2);
    hold on;
    plot(377*Hz,'r','LineWidth',2);
    set(gca,'YLim',ylims); 

fprintf('-------------------------\n')
fprintf('Nx = %d\n',Nx);
fprintf('Lx = %.1f [m]\n',Lx);
fprintf('dx = Lx/Nx = %.1e [m]\n',ds);
fprintf('fs = %. 1e [Hz]\n',fs);
fprintf('dt = ds/fs = %.1e [s]\n',dt);
fprintf('i_lamda = %.1f\n',Nx*c/fs/Lx);
fprintf('lamda   = %.2e [m]\n',c/fs/Lx);
fprintf('max(sigma)*2*pi*f/epsilon_o = %.1e\n',max(sigma)*2*pi*fs/eps0);
fprintf('-------------------------\n')

geometric_ratio = sqrt(mur(ixb-1)*epsr(ixb)/mur(ixb)*epsr(ixb-1));
reflection_coefficient = (1-geometric_ratio)/(1+geometric_ratio);


for iter=1:Niter
    % Source
    if source == 1
        Ey(2) = sin(2*pi*fs*dt*iter);
        Ey_forward(2) = Ey(2);
    end
    Eforward = Ey(2);
    if source == 2
        % Gaussian pulse
        Ey(3) = exp(-((iter-10)/5)^2);
    end
    
    % The next 10 or so lines of code are where we actually integrate Maxwell's
    % equations. All the rest of the program is basically bookkeeping and plotting.
    Hz(1) = Hz(2); % Absorbing boundary conditions for left-propagating waves
    Hz_forward(1) = Hz_forward(2);
    for i=2:Nx-1 % Update H field
      Hz(i) = Hz(i)-am(i)*(Ey(i+1)-Ey(i));
      if i < ixb - 1
        Hz_forward(i) = Hz_forward(i)-am(i)*(Ey_forward(i+1)-Ey_forward(i));
      end
    end
    Ey(Nx) = Ey(Nx-1); % Absorbing boundary conditions for right-propagating waves
    Ey_forward(ixb-1) = Ey_forward(ixb-2); % Absorbing boundary conditions for right-propagating waves

    for i=2:Nx-1 % Update E field
      Ey(i) = as(i)*Ey(i)-ae(i)*(Hz(i)-Hz(i-1));
      if i < ixb -1 
        Ey_forward(i) = Ey_forward(i)-ae(i)*(Hz_forward(i)-Hz_forward(i-1));
      end
    end
    
    Hz_ip(iter) = Hz(ip);
    Ey_ip(iter) = Ey(ip);

    Ey_reflected = Ey(1:(ixb-1))-Ey_forward;
    Ey_pen = Ey(ixb:end);
    
    if (animate || iter == Niter)
        figure(2);hold off;
            plot(Ey(1:(ixb-1)),'b','LineWidth',2);
            
            hold on;
            grid on;
            plot(Ey_forward,'c','LineWidth',2)
            plot(Ey_reflected,'g','LineWidth',2)
            plot(ixb:Nx,Ey_pen,'k','LineWidth',2)
%             plot(377*Hz,'r','LineWidth',2);
%             plot(377*Hz_forward,'m')
            title(sprintf('i_t = %03d; f = %s [Hz]; L_x = %.1f [m]',...
                iter,fstr,Lx));
            if 0  %1
                text(1,ylims(2),sprintf('E_y [V/m] (blue)\n377 H_z [T] (red)'),...
                    'VerticalAlignment','Top');
            elseif animate == 0
%                 legend('E_y [V/m]','377H_z [T]') % Slows down rendering
                legend('E_{y1}','E^{+}_{y1}','E^{+}_{y1}','E_{y2}');
            end
            xlabel('i_x');
            fdtd1d_annotate;
            drawnow;
            %pause(0);
            
    end

    if iter == 1
        fprintf('i_t (Time Step) = %04d',iter);
    end
    if iter > 1
        fprintf('\b\b\b\b');
        fprintf('%04d',iter);
    end

end
fprintf('\n');

% fdtd1d_plot;
