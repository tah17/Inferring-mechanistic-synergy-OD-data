function model = edwards_2drug()

% Symbolic variables
syms x1 x2 x3 x4 x5 x6 x7
syms L b basal d_inv g_p_1 k_1 e_1 g_p_2 k_2 e_2 g_p_12 k_12 e_12

% Parameters
model.sym.p = [L; b; basal; d_inv; g_p_1; k_1; e_1; g_p_2; k_2; e_2; g_p_12; k_12; e_12];

% State variables
model.sym.x = [x1; x2; x3; x4; x5; x6; x7];

% Control vectors (g)
model.sym.g = [];

% Autonomous dynamics (f)
model.sym.xdot = [b*x1*exp(-x1/L); 
    (b/(1+g_p_1))*x2*exp(-x2*(1+e_1)/L) - k_1*x2; 
    k_1*x2; 
    (b/(1+g_p_2))*x4*exp(-x4*(1+e_2)/L) - k_2*x4; 
    k_2*x4;
    (b/(1+g_p_1+g_p_2+g_p_12))*x6*exp(-x6*(1+(e_1+e_2+e_12))/L) - (k_1+k_2+k_12)*x6; 
    (k_1+k_2+k_12)*x6];

% Initial conditions
model.sym.x0 = [250; 250; 0; 250; 0; 250; 0];

% Observables
model.sym.y = [basal; 
    basal + x1*d_inv; 
    basal + (x2 + x3)*d_inv;
    basal + (x4 + x5)*d_inv;
    basal + (x6 + x7)*d_inv];

end