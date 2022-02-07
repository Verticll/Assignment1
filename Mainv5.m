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
    Rx = zeros(1,nElec)-200e-9;
    Ry = 2 * (rand(1, nElec)-0.5);
    Rtheta = (180 * rand(1, nElec))-90;
    RV = Vth * abs(randn(1,nElec));
    %TIME
    t = 0;
    dt = 1e-13;
    TStop = 1000*dt;
    %COUNT
    count = 1;
    nElecCount = 0;
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
    %instantiate variables
    xOld = zeros(1,nElec);
    yOld = zeros(1,nElec);

    
    % randomly place  abunch of particles 1000-10000
    x(1, :) = Rx;
    y(1, :) = Ry * yMax;
    
    % give each particle Vth but with a random direction
    Vx(1:nElec) = RV .* cos(Rtheta);
    Vy(1:nElec) = RV .* sin(Rtheta);
    VTEST = mean(sqrt(Vx.^2 + Vy.^2));
    VTEST1 = mean(RV);
   
    while t < TStop
            

            %Generate normal distribution for scatter tests
            rScat = rand(1,nElecCount);
            scattercount = 0;

            %Scatter test
            for i=1:1:nElecCount
               if rScat(i) <= pScat
                  scattercount = scattercount+1;
                  RV1 = Vth * abs(randn(1, 1));
                  Rtheta1 = 360 * rand(1, 1);
                  Vx(i) = RV1 * cos(Rtheta1);
                  Vy(i) = RV1 * sin(Rtheta1);
               end
            end


            %Get positions
           
            xOld(1:nElecCount) = x(1:nElecCount);
            yOld(1:nElecCount) = y(1:nElecCount);
            x(1:nElecCount) = x(1:nElecCount) + (Vx(1:nElecCount) .* dt);
            y(1:nElecCount) = y(1:nElecCount) + (Vy(1:nElecCount) .* dt);

            %Get plot arrays using new - old
            xPlot = [xOld(:) x(:)];
            yPlot = [yOld(:) y(:)];
            
            %Iterate time
            t  = t + dt;              
            
            %Reflect on the Ymax
            for i=1:1:nElecCount
               if y(i) <= -yMax || y(i) >= yMax
                  Vy(i) = Vy(i) * -1; 
               end
            end
            for i=1:1:nElecCount
                if x(i) <= -xMax || x(i) >= xMax
                  Vx(i) = Vx(i) * -1; 
               end
            end

            %Checking box boundaries
            for i=1:1:nElecCount
                for j=1:1:nElecCount
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
            ScatterAvg = scattercount/nElecCount; % get current chance to scatter
            TempPlot(count,:) = TAvg;  % Vth = sqrt(C.kb * 300/ (0.26*C.m_0));
            scatterPlot(count,:) = ScatterAvg; % store array of times
            timePlot(count,:) = [t]; % store array of times
            
            subplot(2,2,1);
            hold on
            %subplot(2,1,1),plot(x, y, 'bo', 'markers',4,'MarkerFaceColor', 'b');
            subplot(2,2,1),plot(Box(1:5,1),Box(1:5,2),'k');
            if(nElecCount >= 1)
            subplot(2,2,1),plot(xPlot(1,1:2), yPlot(1,1:2),'b');
            end
            if(nElecCount >= 14)
            subplot(2,2,1),plot(xPlot(11,1:2), yPlot(11,1:2),'r');
            end
            if(nElecCount >= 24)
            subplot(2,2,1),plot(xPlot(21,1:2), yPlot(21,1:2),'g');
            end
            if(nElecCount >= 34)
            subplot(2,2,1),plot(xPlot(31,1:2), yPlot(31,1:2),'c');
            end
            if(nElecCount >= 44)
            subplot(2,2,1),plot(xPlot(41,1:2), yPlot(41,1:2),'y');
            end
            if(nElecCount >= 54)
            subplot(2,2,1),plot(xPlot(41,1:2), yPlot(41,1:2),'m');
            end
            
            
            
            
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
            
            if nElecCount < 1000
                nElecCount = nElecCount+10;
            end
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

