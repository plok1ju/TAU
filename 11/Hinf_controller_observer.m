function [Q, K, P, L] = Hinf_controller_observer(A, B_1, B_2, C_1, C_2, D_1, D_2, gamma)
    [Q, K, ~] = icare(A, B_2, C_2'*C_2, D_2'*D_2, [], [], gamma^(-2)*B_1*B_1');
    K = -K;
    [P, ~, ~] = icare(A', C_1', B_1*B_1', D_1*D_1', [], [], gamma^(-2)*C_2'*C_2);
    L = -P*(eye(size(Q)) - gamma^(-2)*Q*P)^-1*(C_1+gamma^(-2)*D_1*B_1'*Q)'*(D_1*D_1')^-1;
end