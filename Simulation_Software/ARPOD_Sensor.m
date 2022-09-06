classdef ARPOD_Sensor
    methods (Static)
        function z_t = measure2D(state)
            x = state(1);
            y = state(2);
            theta = state(3);

            theta_pos = atan(y/x);
            dist = sqrt(x*x+y*y);
            theta_body = theta;
            z_t = [theta_pos; theta_body; dist];
        end
        function jacobian = measureJacobian2D(state)
            x = state(1);
            y = state(2);
            theta = state(3);


            jacobian = zeros(3,6);
            %dArctan(y/x)
            partialX = -y/(x*x+y*y);
            partialY = x/(x*x+y*y);
            jacobian(1,1) = partialX;
            jacobian(1,2) = partialY;

            %sqrt(x^2+y^2)
            norm = sqrt(x*x+y*y);
            partialX = x/norm;
            partialY = y/norm;
            jacobian(3,1) = partialX;
            jacobian(3,2) = partialY;

            %theta_C
            jacobian(2,3) = 1;
        end
        function z_t = sense2D(state,noise_model, phase)
                z_t = ARPOD_Sensor.measure2D(state);
                z_t = z_t + noise_model();
            if phase == 1
                z_t = z_t(1:2,:); % do not include distance
            end
        end
    end
end