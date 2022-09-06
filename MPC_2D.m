classdef MPC_2D
    methods (Static)
        function [H,f] = setupQuadraticCost(Q,R,n_horizon)
            %{
                Parameters:
                ------------
                    Q: 6x6 diagonal matrix that weighs the states in cost
                    function
                    R: 2x2 diagonal matrix that weights the control inputs in
                    cost function
                    n_horizon: scalar horizon that MPC spans over.
                Description:
                ------------
                    Sets up Quadratic Cost function for the MPC as 
                    z^T*H*z + f^T*z
                    Makes H sparse for program to work faster.
                    (works since H is just diagonal)
            %}
            H = [];
            for i = 1:n_horizon
                H = [H, diag(Q).'];
                if i ~= n_horizon
                    H = [H, diag(R).'];
                end
            end
            ndim = length(H);
            f = zeros(ndim,1);
            H = sparse(diag(H));
        end
        function [Aeq, beq] = setupLinearConstraints(traj0, tstep, n_horizon)
            %{
                Parameters:
                ------------
                    traj0: initial trajectory
                    tstep: scalar timestep in seconds
                    n_horizon: scalar horizon that spans MPC
                Description:
                ------------
                    Setups the equality constraints of MPC
                    This is namely the dynamics Ax_k+Bu = x_k+1
                    These matrices are taken from the Discretized HCW
                    equations.
            %}
            %setting up equality constraint (dynamics)
            theta = traj0(3);
            [A,B] = ARPOD_Dynamics.HCW_ZOH_2D(theta, tstep);
            Aeq = [];
            row0 = [];
            for i = 1:n_horizon-1
                row = [];
                for j = 1:n_horizon
                    if j == i
                        row = [row, A,B];
                    elseif j == i+1
                        if i ~= n_horizon-1
                            row = [row, -eye(6), zeros(6,2)];
                        else
                            row = [row, -eye(6)];
                        end
                    elseif j == n_horizon
                        row = [row, zeros(6,6)];
                    else
                        row = [row, zeros(6,6), zeros(6,2)];
                    end
                end

                Aeq = [Aeq;row];
                %setting up the x0 = xbar0 to ensure first term doesn't
                %change
                if i == 1
                    row0 = [row0, eye(6), zeros(6,2)];
                else
                    row0 = [row0, zeros(6,6), zeros(6,2)];
                end
            end
            row0 = [row0, zeros(6,6)];
            Aeq = [Aeq;row0];
            [m,n] = size(Aeq);
            beq = zeros(m,1);
            beq(m-5:m,:) = traj0;
        end
        function [ub,lb] = setupControlInputBoundaries(n_horizon, mass)
            %{
                Parameters:
                ------------
                    n_horizon: scalar horizon that spans MPC
                Description:
                ------------
                    limits the control inputs according to the ARPOD
                    Benchmark paper.
                    ub and lb are specifically made for MATLAB's
                    optimization functions.
            %}
            ubar = ARPOD_Mission.ubar / mass;
            wbar = ARPOD_Mission.wbar;
            ub = [];
            lb = [];
            for i = 1:n_horizon-1
                ub = [ub; inf;inf;inf;inf;inf;inf;ubar;wbar];
            end
            ub = [ub; inf;inf;inf;inf;inf;inf];
            for i = 1:n_horizon-1
                lb = [lb; -inf;-inf;-inf;-inf;-inf;-inf;0;-wbar];
            end
            lb = [lb; -inf;-inf;-inf;-inf;-inf;-inf];
        end
        function LOS = setupLinearLOS(index, n_horizon)
            %{
                Parameters:
                ------------
                    index: index of the trajectory to setup LOS matrix for.
                    n_horizon: scalar horizon that spans MPC
                Description:
                ------------
                    Sets up the LOS constraint for one trajectory in the
                    whole span of trajectories.
                    This is specifically done to easily create LOS
                    constraints for QP and nonlcon
            %}
            theta1 = ARPOD_Mission.theta * pi / 180;
            LOS_mtx = [ sin(theta1/2), cos(theta1/2),0; 
                    sin(theta1/2), -cos(theta1/2),0];
            LOS_mtx = [LOS_mtx, zeros(2,3)];
            LOS = [];
            for i = 1:n_horizon
                if i == index
                    mtx = LOS_mtx;
                    if i ~= n_horizon
                        mtx = [mtx, zeros(2,2)];
                    end
                elseif i ~= n_horizon
                    mtx = zeros(2,8);
                else
                    mtx = zeros(2,6);
                end
                LOS = [LOS, mtx];
            end
        end
        function [xs,us] = extractOptVector(vector, n_horizon)
            %{
                Parameters:
                -----------
                    vector: large vector of a certain form below
                    n_horizon: scalar horizon that spans MPC
                Description:
                ------------
                    vector is setup like so
                    [ x_0
                      u_0
                      x_1
                      u_1
                      ....
                      u_n-1
                      x_n ]
                    This is because of how desirable this format is when
                    performing optimization. 
                    So, this function seeks to split this into just x and
                    u which makes it easier to work with when trying to
                    graph or collect data for benchmarking. 
            %}
            xs = [];
            us = [];
            for i = 1:n_horizon
                x_id = 8*(i-1)+1;
                xs = [xs, vector(x_id:x_id+5)];
                if i ~= n_horizon
                    u_id = 8*(i-1)+1+6;
                    us = [us, vector(u_id:u_id+1)];
                end
            end
        end
        function [xs,us, xstar] = optimizeLinear(traj0, x0, Q, R, n_horizon, tstep, mass, phase)
            %{
                Parameters:
                ------------
                    traj0: 6x1 current trajectory of spacecraft
                    Q: 6x6 weighting matrix for states.
                    R: 2x2 weighting matrix for control inputs.
                    n_horizon: scalar horizon that spans MPC
                    tstep: scalar timestep in seconds.
                    phase: current phase of the spacecraft. 
                Description:
                ------------
                    A combination of all of the algorithms written in this
                    file that sets up an objective function with dynamic
                    constraints, control constraints, and docking
                    constraints. 
                    The endgoal is a fully realized trajectory generator
                    that allows spacecraft to rendezvous and dock
                    successfully.
            %}

            %setting up quadratic cost function J = sum( x^TQx + u^TRu )
            [H,f] = MPC_2D.setupQuadraticCost(Q,R,n_horizon);
            % cost function -> z^THz + f^Tz

            
            if phase == 2
                %f = ChaserMPC.setupLOSRelaxation(1e10,n_horizon);
                traj0 = traj0 + [5*ARPOD_Mission.rho_d/6;0;0;0;0;0];
            end
            

            %setting up equality constraint (dynamics)
            [Aeq, beq] = MPC_2D.setupLinearConstraints(traj0, tstep, n_horizon);

            %setting up bound constraints
            [ub,lb] = MPC_2D.setupControlInputBoundaries(n_horizon, mass);

            % phase 2: only ensure last horizon reaches LOS
            % phase 3: all states are in LOS + velocities are bounded.
            vbar = ARPOD_Mission.Vbar;
            A = [];
            b = [];
            if phase == 3
                for i = 1:n_horizon
                    LOS = MPC_2D.setupLinearLOS(i,n_horizon);
                    A = [A; LOS];
                end
                for i = 1:9:9*n_horizon-3
                    ub(i+3:i+5) = [vbar;vbar;vbar];
                    lb(i+3:i+5) = [-vbar;-vbar;-vbar];
                end
            end
            [n,m] = size(A);
            b = zeros(n,1);

            if x0 == -1
                x0 = zeros(length(H), 1);
            end

            if phase == 3
                options = optimoptions(@quadprog, 'Algorithm', 'active-set');
            else
                options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex');
            end
            final = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options);
            %final = quadprog(H,f,A,b,Aeq,beq,lb,ub);

            [xs,us] = MPC_2D.extractOptVector(final, n_horizon);
            xstar = final;
        end
        function [xstar,u] = controlMPC(traj0, x0, Q, R, n_horizon, tstep, mass, phase)
            %{
                Parameters:
                -----------
                    Same parameters as MPC optimize function
                    With addition of mpc_option
                Description:
                ------------
                    Blackbox controller for MPC (all the opt functions
                    above).
            %}
            [xs,us,xstar] = MPC_2D.optimizeLinear(traj0, x0, Q,R,n_horizon,tstep,mass,phase);
            u = us(1:2,1); %return first control input
        end
    end
end