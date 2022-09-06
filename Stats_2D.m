classdef Stats_2D
    properties
        fuelConsumed = 0;
        currTraj = [0;0;0;0;0;0];
        
        trackTraj = [];
        trackEstTraj = [];
        trackU = [];
        trackFuelConsumption = [];
        trackPhase = [];
        timestamps = [];
        total_steps = 0;

        estimation_ts = [0];
    end
    methods
        function obj = initBenchmark(obj, traj0)
            obj.currTraj = traj0;
            obj.trackFuelConsumption = [0];
            obj.timestamps = [0];
            obj.trackTraj = [obj.trackTraj, traj0];
            obj.trackEstTraj = [obj.trackEstTraj, traj0];
            obj.total_steps = 1;
        end
        function obj = updateBenchmark(obj, u, mass, trueTraj, estTraj, tstep)
            obj.trackU = [obj.trackU, u];

            fuel = obj.trackFuelConsumption(length(obj.trackFuelConsumption));
            % fuel -> newton seconds or kg*m/s or N*s
            fuel = fuel + u(1) * mass * tstep; 
            obj.trackFuelConsumption = [obj.trackFuelConsumption, fuel];

            obj.currTraj = trueTraj;
            obj.trackEstTraj = [obj.trackEstTraj, estTraj];
            obj.trackTraj = [obj.trackTraj, obj.currTraj];
            obj.total_steps = obj.total_steps + 1;
            obj.timestamps = [obj.timestamps, obj.timestamps(length(obj.timestamps)) + tstep];

        end
        function obj = graphLinear(obj, theta1)
            %{
                Graph Constituents
                ------------------
                    - Chaser Start
                    - Chaser End
                    - Chaser Trajectory
                    - Origin 'x' (target position)
                    - 4th rendezvous location 'x'
                    - LOS boundary (linear pyramid)
                    - Phase 2 boundary
                    - Phase 3 boundary
            %}
            set(groot,'defaultAxesTickLabelInterpreter','latex');  
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            
            fsize = 20; % figure fontsize
            lw = 2; % linewidth


            figure(1) % true trajectories and estimates in 2D planes
            set(gcf,'Position', [10 10 800 800])
            plot(obj.trackTraj(1,:), obj.trackTraj(2,:), 'r','LineWidth',lw)
            hold on
            plot(obj.trackEstTraj(1,:), obj.trackEstTraj(2,:), 'b--','LineWidth',lw)
%             title('Chaser Trajectory 2D')
            xlabel("$x$ [km]")
            ylabel("$y$ [km]")
            set(gca,'fontsize',fsize)
            legend('true','estimate','Location','southeast')
            set(gca, 'TickLabelInterpreter', 'latex')
            grid

            return;
        end
        function totalMSE = getError(obj)
            MSEperStep = sum((obj.trackTraj - obj.trackEstTraj).^2);
            totalMSE = sum(MSEperStep);
        end
    end
end