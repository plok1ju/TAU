% система
A = [-9,0,-10;-4,-1,-6;6,-2,5];
B = [0; 0; 0];
C = [1, 0, 1];
D = 0;
P = [-1, -3-1i, -3+1i;-1,-1-1i,-1+1i;1,2,2];
sys = ss(A, B, C, D);
x_1 = [4; 3; -3];
t_1 = 3;
% матрица наблюдения 
V = obsv(sys);

% ранг матрицы наблюдения 
rank_V = rank(V);

% собственные числа А
e = eig(A);
% С в жордановой форме
CP = C*P;

% граммиан
opt = gramOptions('TimeIntervals',[0 t_1]);
GR = gram(sys, 'o', opt);

syms t real;

exp_a = expm(A'*(t_1-t));
exp_a_simplify = simplify(exp_a, "Steps", 100);

u_x_1 = B'*exp_a_simplify/P*x_1;
u_x_1_simplify = simplify(u_x_1);