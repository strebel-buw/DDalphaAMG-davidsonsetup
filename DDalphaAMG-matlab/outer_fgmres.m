function [x,res] = outer_fgmres(A,b,block,settings,setup_flag,x)
% ------------------------------------------------------------
% GMRES
restart = settings.restart;
tol = settings.tol;
maxrestart = settings.maxrestart;

n = length(b);
if nargin < 6 
    x = zeros(n,1);
end
res = norm(b-A*x);
e_1 = zeros(restart+1,1);
e_1(1) = 1;
iter = 0;

% ------------------------------------------------------------

for k = 1:maxrestart
        
    r = b - A*x;
    beta = norm(r);
    
    if maxrestart > 2 && settings.disp == 1
        fprintf(['  sweep ',num2str(k),', ']),
        fprintf(['  current relres: ',num2str(beta/norm(b)),'\n']),
    end
    
    V = zeros(n,restart+1);
    Z = zeros(n,restart);
    V(:,1) = r/beta;
    H = zeros(restart+1,restart);
    
    for j = 1:restart
        % -----------------------------------------------
        % Arnoldi MGS
        % right preconditioning
        iter = iter + 1;
        
        if settings.multigrid == 1 
            Z(:,j) = amg(A,V(:,j),block,settings,setup_flag);
        else
            Z(:,j) = V(:,j);
        end

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
        % -----------------------------------------------
        % Solve least square problem
        % argmin_y || beta e_1 - H y ||
        [~,R] = qr([H(1:j+1,1:j),beta*e_1(1:j+1)]);
        % -----------------------------------------------
        res = [res,abs(R(j+1,j+1))];
        
        if (settings.disp == 1) 
            fprintf('relres after %d iterations: %e\n',iter,res(end)/res(1))
        end
        if res(end)/res(1) < tol || j == restart
            % -----------------------------------------------
            % compute x_0 = x_0 + V*y_m
            y = R(1:j,1:j)\R(1:j,j+1);
            x = x+Z(:,1:j)*y;
            break;
            % -----------------------------------------------
        end
    end
    
    if res(end)/res(1) < tol || H(j+1,j) == 0
        break;
    end 
end
