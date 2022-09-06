% gravitational constant in km^2/s^2 from chance-constr MPC

classdef ARPOD_Dynamics
    %{
        These are dynamics for the underactuated spacecraft

    %}
    properties 
        traj
    end
    methods (Static)
        function [A,B] = HCW_ZOH_2D(theta,T)
            mu_GM = ARPOD_Mission.mu;
            R = ARPOD_Mission.a;

            n = sqrt(mu_GM / (R.^3) );
            A = zeros(6,6);
            B = zeros(6,2);

            S = sin(n * T);
            C = cos(n * T);

            A(1,:) = [4-3*C,0,0,S/n,2*(1-C)/n,0];
            A(2,:) = [6*(S-n*T),1,0,-2*(1-C)/n,(4*S-3*n*T)/n,0];
            A(3,:) = [0,0,1,0,0,T];
            A(4,:) = [3*n*S,0,0,C,2*S,0];
            A(5,:) = [-6*n*(1-C),0,0,-2*S,4*C-3,0];
            A(6,:) = [0,0,0,0,0,1];

            B(1,:) = [(1-C)/(n*n)*cos(theta) + (2*n*T-2*S)/(n*n)*sin(theta), 0];
            B(2,:) = [-(2*n*T-2*S)/(n*n)*cos(theta) + ((4*(1-C)/(n*n))-(3*T*T/2))*sin(theta), 0];
            B(3,:) = [0, 0.5*T*T];
            B(4,:) = [S/n*cos(theta) + 2*(1-C)/n*sin(theta), 0];
            B(5,:) = [-2*(1-C)/n*cos(theta) + (4*S/n - 3*T)*sin(theta), 0];
            B(6,:) = [0,T];
        end
        
        function dtraj = nonlinearMotion2D(traj0, mu_GM, R, u)
            %{
                u = [u_mag, u_wheel]
            %}
            x = traj0(1);
            y = traj0(2);
            theta = traj0(3);

            xdot = traj0(4);
            ydot = traj0(5);
            thetadot = traj0(6);

            u_mag = u(1);
            u_wheel = u(2);


            n = sqrt(mu_GM / (R.^3)); %orbital velocity
            const = ((R+x).^2 + y.^2).^(3/2); 


            xdotdot = 2*n*ydot + n*n*(R+x) - mu_GM*(R+x)/const + u_mag*cos(theta);
            ydotdot = -2*n*xdot + n*n*y - mu_GM*y/const + u_mag*sin(theta);
            thetadotdot = u_wheel;

            dtraj = [xdot;ydot;thetadot;xdotdot;ydotdot;thetadotdot];
        end
        function traj = nonlinearMotionSolver2D(traj0,mu_GM,R,u,tstep)
            motion = @(t,traj) ARPOD_Dynamics.nonlinearMotion2D(traj,mu_GM,R,u);
            if tstep < 1
                [ts,trajs] = ode45(motion, 0:tstep/10:tstep, traj0);
            else
                [ts,trajs] = ode45(motion, 0:tstep, traj0);
            end
            traj = transpose(trajs(end,:));
        end
    end
end