%given data
A = [-9,0,-10;-4,-1,-6;6,-2,5];
B = [0; 0; 0];
C = [1, 0, 1];
D = 0;
sys = ss(A, B, C, D);
check = isstable(sys);
t_1 = 3;

% Observabillity matrix
U = obsv(sys);
rank_U = rank(U);


% observabillity gramian for time t_1
% function gram cannot be use because sym is unstable dynamics
syms t real;
exp_a_1 = simplify(expm((A')*t), "Steps", 100);
exp_a_2 = simplify(expm(A*t), "Steps", 100);
f_g = exp_a_1*(C')*C*exp_a_2;
Gr_t_1 = int(f_g, t, 0, t_1);
Gr_t_1 = double(Gr_t_1);
eig_Gr = eig(Gr_t_1);
 
% initial conditions to lead state vector along the trajectory y(t) during time t_1
y_t = exp(-3*t)*(-3*cos(2*t) - 2*sin(2*t)); 
exp_a = expm(A'*t);
exp_a_simplify = simplify(exp_a, "Steps", 100);
f = exp_a_simplify*C'*y_t;
integral_f = int(f, t, 0, t_1);
x_0 = double((inv(Gr_t_1)*integral_f));