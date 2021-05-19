function [x,res] = fgmres(A,b,restart,tol,maxrestart,disp,M,x)
% ------------------------------------------------------------
% GMRES

if isempty(restart)
    restart = maxrestart;
    maxrestart = 1;
end

n = length(b);
if nargin<8
    x = zeros(n,1);
end

if ( isa(A,'function_handle') )
    res = norm(b-A(x));
else
    res = norm(b-A*x);
end

e_1 = zeros(restart+1,1);
e_1(1) = 1;

% ------------------------------------------------------------

for k = 1:maxrestart
    
    if ( isa(A,'function_handle') )
        r = b - A(x);
    else
        r = b - A*x;
    end
    
    beta = norm(r);
    
    if maxrestart > 2 && disp == 1
        fprintf(['sweep ',num2str(k),', ']),
        fprintf(['current relres: ',num2str(beta/norm(b)),'\n']),
    end
    
    V = zeros(n,restart+1);
    Z = zeros(n,restart);
    V(:,1) = r/beta;
    H = zeros(restart+1,restart);
    
    for j = 1:restart
        % -----------------------------------------------
        % Arnoldi MGS
        % right preconditioning
        if nargin >= 7 && isa(M,'function_handle')
            Z(:,j) = M(V(:,j));
        else
            Z(:,j) = V(:,j);
        end
        if ( isa(A,'function_handle') )
            w = A(Z(:,j));
        else
            w = A*Z(:,j);
        end
        for i=1:j
            H(i,j) = V(:,i)'*w;
            w = w - H(i,j)*V(:,i);
        end
        H(j+1,j) = norm(w);
        if H(j+1,j) ~= 0
            V(:,j+1) = w/H(j+1,j);
        end
        % -----------------------------------------------
        % Solve least square problem
        % argmin_y || beta e_1 - H y ||
        [Q,R] = qr([H(1:j+1,1:j),beta*e_1(1:j+1)]);
        % -----------------------------------------------
        res = [res,abs(R(j+1,j+1))];
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
