function [x,block] = kcycle(A,r,block,settings,setup_flag)

n = length(r);
cycle = settings.cycle;

e_1 = zeros(cycle+1,1);
e_1(1) = 1;
beta = norm(r);
    
V = zeros(n,cycle+1);
Z = zeros(n,cycle);
V(:,1) = r/beta;
H = zeros(cycle+1,cycle);
    
for j = 1:cycle
    % -----------------------------------------------
    % Arnoldi MGS
    % right preconditioning
    [Z(:,j),block] = amg(A,V(:,j),block,settings,setup_flag);
    w = A*Z(:,j);
    for i=1:j
        H(i,j) = V(:,i)'*w;
        w = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    if abs(H(j+1,j)) > 1e-15
        V(:,j+1) = w/H(j+1,j);
    else
        break
    end
end
% -----------------------------------------------
% Solve least square problem
% argmin_y || beta e_1 - H y ||
% compute x = Z*y_m
[~,R] = qr([H(1:j+1,1:j),beta*e_1(1:j+1)]);
x = Z(:,1:j)*(R(1:j,1:j)\R(1:j,j+1));
