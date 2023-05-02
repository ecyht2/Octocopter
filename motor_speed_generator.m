function [w] = motor_speed_generator(input,k_T,k_M,L)

Ainv = [
    1/(8*k_T), 1/(2*k_T), 0, -1/(8*k_M);
    L/(4*k_T*(L + 1)), -1/(k_T*(2^(1/2)*L + 2^(1/2))), -1/(k_T*(2^(1/2)*L + 2^(1/2))), L/(4*k_M*(L + 1));
    1/(8*k_T), -1/(2*k_T*L), 0, -1/(8*k_M);
    1/(4*k_T*(L + 1)),  1/(k_T*(2^(1/2)*L + 2^(1/2))),  1/(k_T*(2^(1/2)*L + 2^(1/2))), 1/(4*k_M*(L + 1))
    ];

%input = [Uz; U_phi; U_theta; U_psi];

w = Ainv*input;

end