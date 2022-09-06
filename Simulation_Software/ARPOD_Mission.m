classdef ARPOD_Mission
    properties (Constant)
        t_e = 14400; % eclipse time (in seconds)
        t_f = 43200; % total mission duration (in seconds)
        rho_r = 1; % maximum distance for range measurements (1 km)
        rho_d = 0.1; % docking phase initial radius (0.1 km)
        m_t = 2000; % mass of target (2000 kg)
        m_c = 500; % mass of chaser (500 kg)
        mu = 398600.4; %earth's gravitational constant (in km^3/s^2)
        a = 42164; % semi-major axis of GEO (42164 km)
        Vbar = 5 * 10.^(-5); % max closing velocity while docking (in km/s)
        theta = 60; % LOS Constraint angle (in degrees)
        c = [-1;0;0]; % LOS cone direction
        x_docked = [0;0;0;0;0;0]; % docked position in km, and km/s
        x_relocation = [0;20;0;0;0;0]; %relocation position in km, km/s
        x_partner = [0;30;0;0;0;0]; %partner position in km, km/s
        ubar = 10; % max 10 Newtons

        wbar = 1; % max 1 Nkm
    end
    methods (Static)
        function inLOS = isInsideLOS_2D(traj)
            theta1 = ARPOD_Mission.theta * pi / 180; %docking angle in radians
            LOS_mtx = [ sin(theta1/2), cos(theta1/2); 
                    sin(theta1/2), -cos(theta1/2)];
            xy = [traj(1);traj(2)];
            b = LOS_mtx*xy;
            if b(1) <= 0 && b(2) <= 0
                inLOS = 1;
            else
                inLOS = 0;
            end
        end
        
        function phase = calculatePhase2D(traj)
            norm = traj(1:2,:);
            norm = sqrt(sum(norm.^2));
            if (norm > ARPOD_Mission.rho_r)
                % ARPOD phase 1: Rendezvous w/out range
                phase = 1;
            elseif ( (norm > ARPOD_Mission.rho_d) || (ARPOD_Mission.isInsideLOS_2D(traj) == 0) ) 
                % ARPOD phase 2: Rendezvous with range
                phase = 2;
            else 
                %ARPOD phase 3: Docking
                phase = 3;
            end
        end

    end
    properties
        traj % current trajectory of mission
        sensor % current measurements from sensors
        phase % current phase of mission
        inLOS % boolean as whether it is in or not

        los_type %type of los, linear or nonlinear
    end
    methods
        function obj = initMission2D(obj, traj)
            obj.traj = traj;
            obj.phase = ARPOD_Mission.calculatePhase2D(traj);
            obj.inLOS = ARPOD_Mission.isInsideLOS_2D(traj);
        end
        function obj = nextStep2D(obj, control, system_noise, sensor_noise, tstep)
            obj.traj = ARPOD_Dynamics.nonlinearMotionSolver2D(obj.traj, ARPOD_Mission.mu, ARPOD_Mission.a, control, tstep);
            obj.traj = obj.traj + system_noise();

            obj.phase = ARPOD_Mission.calculatePhase2D(obj.traj);
            obj.sensor = ARPOD_Sensor.sense2D(obj.traj, sensor_noise, obj.phase);
            obj.inLOS = ARPOD_Mission.isInsideLOS_2D(obj.traj);
        end
    end
end
