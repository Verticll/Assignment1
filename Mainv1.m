% Question 1 Thermal Velocity at T=300K?
% Vth = sqrt(C.kb * 300K / C.m_0)

% Question 2 mean free path?

% Write a program as follows
% at interval of dt update locations using newtons laws of motion. 
% time step should be 1/100th of the region size
% simulate for 1000 steps
% trace some of the trajectories with plot
% show 2d plot of all or some particles
% for the y boundarie reflect particles 
% for the x boundarie jump to opposite side
% show temperature on plot
% use arrays for position and velocity

function [] = Mainv1(nElec)
global C

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      %metres (32.1740 ft) per s²
    
    %CONSTANTS
    m_Si = 4.6637066e-23;
    m_elec = 0.26 * C.m_0;
    Vth = sqrt(C.kb * 300/ m_Si);
    %RANDOM VALUES
    Rx = 2 * (rand(1, nElec)-0.5);
    Ry = 2 * (rand(1, nElec)-0.5);
    Rtheta = 360 * rand(1, nElec);
    %TIME
    t = 0;
    dt = 100e-11;
    TStop = 100e-8;
    %COUNT
    count = 1;
    %BOUNDARIES
    xMax = 200e-9;
    yMax = 100e-9;
    Limits = [-xMax +xMax -yMax +yMax];
    LimitsTime = [0 TStop 0 12];
    %FLAGS
    xFlag = 0;
    
    % randomly place  abunch of particles 1000-10000
    x(1, :) = Rx * xMax;
    y(1, :) = Ry * yMax;
    
    % give each particle Vth but with a random direction
    Vx(1:nElec) = Vth * cos(Rtheta);
    Vy(1:nElec) = Vth * sin(Rtheta);
    meanVy = mean(Vy);
    meanVx = mean(Vx);
   
    while t < TStop
            %Get positions
            xOld = x;
            yOld = y;
            x = x + (Vx .* dt);
            y = y + (Vy .* dt);
            %Iterate time
            t  = t + dt;  
            
            %Get trajectories
            xCaptured(count,:) = [x(1:5)];
            yCaptured(count,:) = [y(1:5)];
            
            %Reflect on the Ymax
            for i=1:1:nElec
               if y(i) <= -yMax || y(i) >= yMax
                  Vy(i) = Vy(i) * -1; 
               end
            end
            %translation on the Xmax
            for i=1:1:nElec
               if x(i) <= -xMax
                  x(i) = x(i) + 2 * xMax;
                  
               else if x(i) >= xMax
                  x(i) = x(i) + 2 * -xMax;      
                   end
               end
            end
            
            %Temperature calc
            Vavg = sqrt(mean(abs(Vx))^2 + mean(abs(Vy))^2); % get current average Velocity
            VavgPlot(count,:) = [Vavg]; % store array of times and average velocities
            timePlot(count,:) = [t];
            
            subplot(2,1,1),plot(xCaptured, yCaptured);
            hold on
            %subplot(2,1,1),plot(x, y, 'bo', 'markers',4,'MarkerFaceColor', 'b');
            subplot(2,1,1),plot(xCaptured, yCaptured);
            %quiver(x,y,Vx,Vy);
            hold off
            axis(Limits);
            xlabel('x');
            ylabel('y');
            grid on
            
            subplot(2,1,2); %plot Avg Velocity
            plot(t,Vavg,'v','linewidth', 2);
            hold on
            subplot(2,1,2), plot(timePlot,VavgPlot,'linewidth', 2);
            hold off
            axis(LimitsTime);
            xlabel('Time');
            ylabel('Velocity');
            grid on
            
            count = count + 1;

            pause(0.0001)
            
    end
end



% Assign random velocity to start
% Use Maxwell-Blotzmann distribution for each velocity component
% ensure avg speed is vth
% plot the distribution in a histogram
% model the scattering of the electrons using an exponential scattering
% probability
% When they scatter rethermalize them with the MB distribution
% What happens to the avg temp 
% Measure mean free path and meant time between collisions

% add a bottle neck that reflects
% make all boundaries capable of specular of diffusive
% electron density map
% temp map with colors

% curved surface implementation

% injection

