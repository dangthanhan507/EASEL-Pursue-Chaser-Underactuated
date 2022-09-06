close all
clear
clc
rng(1);

traj = [-10;-10;pi/2;-0.01;0.01;0];
total_time = 1000;
tstep = 0.33; % update every second
phase = ARPOD_Mission.calculatePhase2D(traj);

mpc_horizon = 5;
scale_mpcQ = 1;
scale_mpcR = 1;
mpc_Q = scale_mpcQ*[1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        0,0,0,100,0,0;
        0,0,0,0,100,0;
        0,0,0,0,0,0];
mpc_R = scale_mpcR*diag([1,0]);


%EKF
stateEstimator = EKF_2D;
stateEstimator = stateEstimator.initEKF(traj, 1e-40*eye(6)); %really trust initial estimate.

process_noise = 0;

% tunable parameters
if process_noise == 0
    seQ = 1e-20*diag([1,1,1,1,1,1]);
    seR = diag([0.001,0.0001,0.01]);
else
    seQ = diag(zeros(1,6)+process_noise);
    seR = diag([0.001,0.1,0.01]);
end

mission = ARPOD_Mission;
mission = mission.initMission2D(traj);
estTraj = traj;

stats = Stats_2D;
stats = stats.initBenchmark(traj);

for i = tstep:tstep:total_time
    disp(i/tstep)

    noise_noise = 0;
    process_noise_noise = 0;
    if phase == 1
        noiseQ = @() mvnrnd([0;0;0;0;0;0], [0,0,0,0,0,0] + process_noise).';
        noiseR = @() mvnrnd([0;0;0], [0.001, 0.0001, 0.01] + noise_noise).';
    elseif phase == 2
        noiseQ = @() mvnrnd([0;0;0;0;0;0], [0,0,0,0,0,0] + (process_noise + process_noise_noise)*0.1).';
        noiseR = @() mvnrnd([0;0;0], [0.001, 0.0001, 0.01] + noise_noise).';
    else
        noiseQ = @() mvnrnd([0;0;0;0;0;0], [0,0,0,0,0,0] + (process_noise + process_noise_noise)*0.01).';
        noiseR = @() mvnrnd([0;0;0], [0.001, 0.0001, 1e-5] + noise_noise).';
    end
    
    if (i == tstep)
        x0 = -1;
        [x0,u] = MPC_2D.controlMPC(estTraj,x0,mpc_Q,mpc_R,mpc_horizon,tstep,ARPOD_Mission.m_c,phase);
    else
        [A,B] = ARPOD_Dynamics.HCW_ZOH_2D(estTraj(4),tstep);
        x0 = [ x0(9:end,:); [0;0;0]; A*x0(end-5:end)];
        [x0,u] = MPC_2D.controlMPC(estTraj,x0,mpc_Q,mpc_R,mpc_horizon,tstep,ARPOD_Mission.m_c,phase);
    end
    u(2) = sin( -(estTraj(3) - atan(estTraj(2)/estTraj(1))) );
    u(2) = u(2)/tstep/tstep/500;
    
    mission = mission.nextStep2D(u,noiseQ, noiseR, tstep);

    traj = mission.traj;
    phase = ARPOD_Mission.calculatePhase2D(estTraj);
    meas = ARPOD_Sensor.sense2D(estTraj, noiseR, phase);

    stateEstimator = stateEstimator.estimate(u,meas,seQ,seR,tstep,phase);
    estTraj = stateEstimator.state;
    stats = stats.updateBenchmark(u, ARPOD_Mission.m_c, traj, estTraj, tstep);

    if sqrt(estTraj(1).^2 + estTraj(2).^2) < 1e-3
        disp("Docked Early")
        break
    end
end
stats.graphLinear(0);