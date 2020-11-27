Nx = 550;
x = [1:Nx];
lambda = 100;
Mx = 300;  % I am going to keep the boundary at 300

% What you've done is a little more advanced than I was looking for.
% You have computed the fully analytic solution a wave reflected off of a
% boundary at ix = 300. To fully answer the question, you would need to
% do a simulation with the the reflection at ix = 550. The idea of this
% problem, which perhaps was not clear, is to "see" the wave that
% is associated with the VSWR by plotting two waves traveling in opposing 
% directions. 

% I am preallocating my reflected wave to match the domain
E_reflected = NaN(1, Mx);
E_transmit = NaN(1,Nx-Mx);
rho = -0.5;  % I am setting rho here to make it easier to test

figure(1);clf
for i = 1:900
    if i >= Mx
        j = max(1,(2*Mx-i-1)); % maintain the domain
        E_reflected(j:end-1) = E_reflected((j+1):end);
        E_reflected(end) = rho*Ei(Mx);  % Get the value at the boundary 
    end
    
    Ei = cos(2*pi*(x-i)/lambda);
    % Set values to right of wave front to NaN so they won't be plotted.
    Ei(i+1:end) = NaN; 

    E_total(1:Mx) = Ei(1:Mx) + E_reflected;
    E_total(isnan(E_total)) = Ei(isnan(E_total));

    if i > Mx
        k = min(i,Nx)-Mx;
        E_transmit(2:k) = E_transmit(1:k-1);
        % based on the fdtd code you gave us before, the frequency changes
        % when it hits a boundary, but I am not able to reproduce that
        E_transmit(1) = E_total(Mx);        
    end
    
    % Plot current time step as light grey.
    if i > Mx
        plot((j+2):Mx,E_total((j+2):end),'k','LineWidth',1,'Color',[1,1,1,0.4]/2);
    end
    % Keep past time steps
    hold on;
    grid on;
    if i > 1
        % Delete previous current time step thick black line
        delete(h)
        delete(h1)
        delete(h2)
        delete(h3)
        delete(p)
    end
    % Plot current time step as thick black line
    h = plot(Ei,'k','LineWidth',2);
    
    h1 = plot(E_reflected, 'r', 'LineWidth',2);
    h2 = plot(E_total, 'b:');
    h3 = plot((Mx+1):Nx,E_transmit, 'g', 'LineWidth',2);
    
    p = patch([Mx, Nx, Nx, Mx],[2, 2, -2, -2],'k');
    set(p,'FaceAlpha',0.1,'EdgeColor','none');
    
    set(gca,'Ylim',[-2,2]);
    set(gca,'Xlim',[1,Nx]);
    %legend('V^{+}');
    % Uncomment the following to hide past time steps
    
    % It is bad form to hard code a solution, but so be it
    
    hold off;
%     if mod(i,300) == 0
%         % Allow early termination of animation
%         input('Continue?');
%     end
    
    drawnow;
end