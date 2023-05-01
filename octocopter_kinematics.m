function [Xdot] = octocopter_kinematics(X, wmotor, m, l, k_T, k_M, lxx, lyy, lzz)
    % Calculates the kinematics of an Octocopter in plus configuration.
    assert(not(isempty(X)) || ~isempty(wmotor));
    assert(size(wmotor, 1) == 8);
    
    % Define constants
    g = 9.81;
    
    J = [lxx, 0, 0; 0, lyy, 0; 0, 0, lzz];
    Jinv = 1 \ J;
    
    % Unpack state variables
    X = num2cell(X);
    [~, ~, ~, u, v, w, phi, theta, psi, p, q, r] = X{:};
    wmotor = num2cell(wmotor);
    [w1, w2, w3, w4, w5, w6, w7, w8] = wmotor{:};
    
    pqr = [p; q; r];
    uvw = [u; v; w];
    
    % Define trigonometric function
    c_phi = cos(phi); s_phi = sin(phi);
    c_theta = cos(theta); s_theta = sin(theta);
    c_psi = cos(psi); s_psi = sin(psi);
    
    % Transformation matrix
    R = [
        c_theta * c_psi, s_phi * s_theta * c_psi - c_phi * s_psi, c_phi * s_theta * c_psi + s_phi * s_psi;
        c_theta * s_psi, s_phi * s_theta * s_psi + c_phi * c_psi, c_phi * s_theta * s_psi - s_phi * c_psi;
        -s_theta, s_phi * c_theta, c_phi * c_theta;
        ];
    
    % Calculate total thrust and torques
    thrust = k_T * (w1^2 + w2^2 + w3^2 + w4^2 + w5^2 + w6^2 + w7^2 + w8^2);
    M_phi = k_T * l * (-w1^2 - w2^2 + w4^2 + w5^2 + w6^2 - w8^2);
    M_theta = k_T * l * (-w2^2 - w3^2 - w4^2 + w6^2 + w7^2 + w8^2);
    M_psi = k_M * (w1^2 - w2^2 + w3^2 - w4^2 + w5^2 - w6^2 + w7^2 - w8^2);

    % Calculate linear accelerations
    pne_dot = R * uvw;
    uvw_dot = cross(pqr, uvw) + [0; 0; thrust / m] - R' * [0; 0; g];

    % Calculate Euler angles rates
    Reuler = [
        1, s_phi * tan(theta), c_phi * tan(theta);
        0, c_phi, -s_phi;
        0, s_phi / c_theta, c_phi / c_theta;
        ];
    euler_dot = Reuler * pqr;

    pqr_dot = Jinv * cross(pqr, J * pqr) + Jinv * [M_phi; M_theta; M_psi];

    % Pack derivative of state variables
    Xdot = [pne_dot; uvw_dot; euler_dot; pqr_dot];
end
