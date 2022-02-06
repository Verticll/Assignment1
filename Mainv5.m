function [] = Mainv5(nElec)
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
    C.e = 2.7182818;                    %wulers number
    
    %CONSTANTS
    Vth = sqrt(C.kb * 300/ (0.26*C.m_0));
    %RANDOM VALUES
    Rx = 2 * (rand(1, nElec)-0.5);
    Ry = 2 * (rand(1, nElec)-0.5);
    Rtheta = 360 * rand(1, nElec);
    RV = Vth * abs(randn(1,nElec));
    %TIME
    t = 0;
    dt = 1e-13;
    TStop = 1000*dt;
    %COUNT
    count = 1;
    %BOUNDARIES
    xMax = 200e-9;
    yMax = 100e-9;
    Limits = [-xMax +xMax -yMax +yMax];
    LimitsTime = [0 TStop 0 300];
    LimitsScatter = [0 TStop 0 1];
    %Scattering probability
    pScat = 1 - C.e^(-(dt/(0.2e-12)));
    scattercount = 0;
    %Drawing of box
    Box = [-40e-9 yMax/2; 40e-9 yMax/2; 40e-9 -yMax; -40e-9 -yMax;-40e-9 yMax/2;];

    
    % randomly place  abunch of particles 1000-10000
    x(1, :) = Rx * xMax;
    y(1, :) = Ry * yMax;
    
    % give each particle Vth but with a random direction
    Vx(1:nElec) = RV .* cos(Rtheta);
    Vy(1:nElec) = RV .* sin(Rtheta);
    VTEST = mean(sqrt(Vx.^2 + Vy.^2));
    VTEST1 = mean(RV);
   
    while t < TStop
            %Generate normal distribution for scatter tests
            rScat = rand(1,nElec);
            scattercount = 0;

            %Scatter test
            for i=1:1:nElec
               if rScat(i) <= pScat
                  scattercount = scattercount+1;
                  RV1 = Vth * abs(randn(1, 1));
                  Rtheta1 = 360 * rand(1, 1);
                  Vx(i) = RV1 * cos(Rtheta1);
                  Vy(i) = RV1 * sin(Rtheta1);
               end
            end


            %Get positions
            xOld = x;
            yOld = y;
            x = x + (Vx .* dt);
            y = y + (Vy .* dt);

            %Get plot arrays using new - old
            xPlot = [xOld(:) x(:)];
            yPlot = [yOld(:) y(:)];
            
            %Iterate time
            t  = t + dt;              
            
            %Reflect on the Ymax
            for i=1:1:nElec
               if y(i) <= -yMax || y(i) >= yMax
                  Vy(i) = Vy(i) * -1; 
               end
            end
            for i=1:1:nElec
               if x(i) <= -xMax
                  x(i) = x(i) + 2 * xMax;
                  xOld(i) = xOld(i) + 2 * xMax;
               else if x(i) >= xMax
                  x(i) = x(i) + 2 * -xMax;
                  xOld(i) = xOld(i) + 2 * xMax;
                   end
               end
            end

            %Checking box boundaries
            for i=1:1:nElec
                for j=1:1:nElec
                    if y(i) <= yMax/2 && -40e-9 <= x(j) && x(j) <= 40e-9
                        if y(i) <= yMax/2
                            Vy(i) = Vy(i) * -1;
                        end
                        if -40e-9 <= x(j) && x(j) <= 40e-9
                            Vx(j) = Vx(j) * -1;
                        end
                    end
                end
            end
            
            %Temperature calc
            VAvg = sqrt(mean(abs(Vx))^2 + mean(abs(Vy))^2); % get current average Velocity
            TAvg = (VAvg.^2 .* (0.26*C.m_0))/C.kb; % get current average Temp in K
            ScatterAvg = scattercount/nElec; % get current chance to scatter
            TempPlot(count,:) = TAvg;  % Vth = sqrt(C.kb * 300/ (0.26*C.m_0));
            scatterPlot(count,:) = ScatterAvg; % store array of times
            timePlot(count,:) = [t]; % store array of times
            
            subplot(2,2,1);
            hold on
            %subplot(2,1,1),plot(x, y, 'bo', 'markers',4,'MarkerFaceColor', 'b');
            subplot(2,2,1),plot(Box(1:5,1),Box(1:5,2),'k');
            subplot(2,2,1),plot(xPlot(1,1:2), yPlot(1,1:2),'b', ...
            xPlot(2,1:2), yPlot(2,1:2),'r',...
            xPlot(3,1:2), yPlot(3,1:2),'g',...
            xPlot(4,1:2), yPlot(4,1:2),'c',...
            xPlot(5,1:2), yPlot(5,1:2),'y',...
            xPlot(6,1:2), yPlot(6,1:2),'m');
            %quiver(x,y,Vx,Vy);
            hold off
            axis(Limits);
            xlabel('x');
            ylabel('y');
            grid on
            
            subplot(2,2,2); %plot Avg Velocity
            plot(t,TAvg,'v','linewidth', 2);
            hold on
            subplot(2,2,2), plot(timePlot,TempPlot,'linewidth', 2);
            hold off
            axis(LimitsTime);
            xlabel('Time');
            ylabel('Temperature in K');
            grid on
            
            subplot(2,2,3); %plot Avg Velocity
            plot(t,ScatterAvg,'v','linewidth', 2);
            hold on
            subplot(2,2,3), plot(timePlot,scatterPlot,'linewidth', 2);
            hold off
            axis(LimitsScatter);
            xlabel('Time');
            ylabel('scattering probability');
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

