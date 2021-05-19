
%%
% Setup driver
function [ block ] = adaptive_setup( A, settings )
    
    % initial setup
    if (isequal(settings.blocked, 'davidson') || isequal(settings.blocked,'hybrid'))
        [block,settings] = init_davidson_setup( A, [], settings );
    else
        [block,settings] = init_setup( A, settings );
    end
    
    fprintf('initial setup part done\n');
    
    % iterative part
    if isequal(settings.blocked, 'davidson')
        block = iterative_davidson_setup(A, block, settings);
    else
        if isequal(settings.type, 'eigs')
            block = eigs_setup(A, block, settings);
        else
            block = iterative_setup(A, block, settings);
        end
    end
    
end


%%
% Eigs-setup for comparison
function block = eigs_setup(A, block, settings)
    
    [block.V,~] = eigs(A'*A,settings.Nvec,'sr');
    block = get_P_Ac( A, block, settings );
    if block.level > 1
        block.coarser_block = eigs_setup( block.coarser_A, block.coarser_block, block.coarser_settings );
    end
end


%%
% initial Davidson setup
function [ block, settings ] = init_davidson_setup( A, V, settings )
    
    settings.coarsegrid = 1;
    block = get_blocks(A,settings);
    settings.coarsegrid = 0;
    
    block.g5 = kron(speye(size(A,1)/settings.dim(1)),diag([-ones(1,settings.dim(1)/2) ones(1,settings.dim(1)/2)]));
    block.kk = 1;
    
    if isempty(V)
        V = ones(size(A,1),1);
    end
    
    block.V = []; block.AV = [];
    block.N = []; block.M = [];
    block.T = []; block.X = [];
    block.lambda = [];
    block = extend_subspace( A, V, block, settings );
    
    for i = 1:settings.Ninit
        
        if (i==1) 
            [block, settings] = davidson_loop( A, block, settings, (i-1)*settings.Nvec, settings.Nvec-1 ); % one step already done
        else
            [block, settings] = davidson_loop( A, block, settings, (i-1)*settings.Nvec, settings.Nvec );
        end
        
        block = clean_up_subspace( block, settings );

        fprintf('level %d -- performed %d non-blocked setup iterations\n',settings.level,i);
        if settings.done == 1 
            fprintf('level %d -- all test vectors reached setup top %e\n',settings.level,settings.setup_tol);
            break
        end
       
    end
    
    %normalize V (for debugging)
    for i=1:settings.Nvec 
        if real(block.V(1,i))<0 
            block.V(:,i) = block.V(:,i)*-1; 
        end 
    end
    
    settings.coarsegrid = 1;
    
    block = get_P_Ac( A, block, settings );
    
    if block.level > 1
        block = get_coarser_settings( block, settings );
        
        PV = [];
        if settings.topdown
            PV = block.P'*block.V;
        end
        [block.coarser_block,block.coarser_settings] = init_davidson_setup( block.coarser_A, PV, block.coarser_settings );
        
    elseif settings.bottomup
        [block.coarser_block.PV,~] = eigs(block.coarser_A' * block.coarser_A, settings.Nbottomup, 'sm');
    end
    
    if settings.bottomup
        block.PV = block.P*block.coarser_block.PV;
        block = extend_subspace( A, block.PV, block, settings );
        block = clean_up_subspace( block, settings );
    end
    
    %ev_visualize( block, settings );
end


%%
% iterative Davidson setup on the finest level
function block = iterative_davidson_setup(A, block, settings)
    settings.coarsegrid = 0;
    if settings.level == settings.Nlevel
        n = settings.Nsetup*settings.Nvec;
    else
        n = settings.Nblocked;
    end
    
    
    for i = 1:settings.Nblocked:n
        
        block = davidson_loop( A, block, settings, i-1, settings.Nblocked );
        
        block = extended_hierarchy_update( A, [], block, settings );
        
        if block.level > 1
            block.coarser_block = iterative_davidson_setup( block.coarser_A, block.coarser_block, block.coarser_settings );
        elseif settings.bottomup
            [block.coarser_block.PV,~] = eigs(block.coarser_A' * block.coarser_A, settings.Nbottomup, 'sm');
        end
        
        if settings.bottomup
            block.PV = block.P*block.coarser_block.PV;
            block = extend_subspace( A, block.PV, block, settings );
            block = clean_up_subspace( block, settings );
        end
        
        ev_visualize( block, settings );
%         block.level
%         [~,I] = sort(abs(block.T));
%         block.T(I(1:10))
%         for k=1:block.kk-1
%             norm(block.g5*A*block.X(:,k) - block.X(:,k)*block.lambda(k))
%         end
        fprintf('level %d -- performed %d davidson iterations\n',settings.level,n);
        if settings.done ==1 
            fprintf('level %d -- all test vectors reached setup top %e\n',settings.level,settings.setup_tol);
            break
        end        
    end

end


%%
function ev_visualize( block, settings )
    if settings.level == settings.Nlevel && 0
        figure(12345+settings.level); clf;
        Ac = block.coarser_block.g5*block.coarser_A;
        plot( sort(real(eig(full(Ac)))), zeros(size(Ac,1),1), 'ob' )
        hold on
        plot( real(diag(block.V'*block.g5*block.AV)), zeros(size(block.M,1),1), 'xr' )
        plot( real(eig(block.N,block.M)), zeros(size(block.M,1),1), 'dk' )
%         plot( real(settings.D), zeros(size(settings.D,1),1), 'g.' )
        hold off
        legend('A_c','V* g5 A V','sigma')
        axis([-0.5 0.5 -0.2 0.2])
        drawnow
%         keyboard
    end
end


%%
function [block,settings] = davidson_loop( A, block, settings, iter, n )
    for j = 1:n
        if settings.simpleres == 0
            [r,block,settings] = get_residual( block, settings, iter+j );
        else            
            r=block.V(:,end);
        end
        if(settings.onD) 
            V = amg(A,r,block,settings,0);
        else
            V = amg(A,block.g5*r,block,settings,0);
        end
        [j, norm(r), norm(V)]
        block = extend_subspace( A, V, block, settings );
        
    end
end


%%
function [r, block, settings] = get_residual( block, settings, i )
    if(settings.onD)
        [Z,T] = eig(block.N);
    else
        [Z,T] = eig(block.N,block.M);
    end
    sigma = diag(T);
    [block.T,I] = sort(abs(sigma));
    %[~,I] = sort(abs(sigma));
    [zz,~] = qr(Z(:,I));
%     zz = Z(:,I);
    for k = 1:size(block.N,1)
        found = 0;
        for j = 1:block.kk-1
            if ( ~isempty(block.lambda) && (abs(sigma(I(k))-block.lambda(j))) < settings.setup_tol )
                found = 1;
                break
            end
        end
        if ~found
            break
        end
    end   
%     sigma(I(k))
    
    W = block.V * zz(:,k);
    AW = block.AV * zz(:,k);
    if(settings.onD)
        r = AW - W*sigma(I(k));
    else
        r = (block.g5*AW) - W*sigma(I(k));
    end
%     norm_r = norm(r)
    if norm(r) < settings.setup_tol 
        fprintf('  lvl%d, iter%d: vec %d/%d (value: %+f%+fi) converged to %e ev-res\n', ...
            settings.level,i,block.kk,size(block.M,1),real(sigma(I(k))),imag(sigma(I(k))), norm(r)/norm(W) )
        block.lambda(block.kk) = sigma(I(k));
        block.X(:,block.kk) = W;
        block.kk = block.kk+1;
        if block.kk > settings.Nvec
            settings.done = 1;
            return
        end
        for k = k:size(block.N,1)
            found = 0;
            for j = 1:block.kk-1
                if ( ~isempty(block.lambda) && (abs(sigma(I(k))-block.lambda(j))) < settings.setup_tol )
                    found = 1;
                    break
                end
            end
            if ~found
                break
            end
        end
        
        W = block.V * zz(:,k);
        AW = block.AV * zz(:,k);
        if(settings.onD)
            r = AW - W*sigma(I(block.kk));
        else
            r = (block.g5*AW) - W*sigma(I(block.kk));
        end
    end
    
end


%%
function block = extend_subspace( A, V, block, settings )
    if ~isempty(V)
        if(settings.onD)
            if size(block.V,1) > 0
                V = V - block.V * (block.V' * V);
            end
            V = V/norm(V);
%             [V,~] = qr(V,0);
            AV = A*V;

            v = [];
            W = AV;
            if size(block.V,1) > 0
                v = block.V' * W;
            end
            y = V' * W;
            block.N = [ block.N, v ; v', y ];            
        else
            if size(block.V,1) > 0
                V = V - block.V * (block.V' * V);
            end
            V = V/norm(V);
%             [V,~] = qr(V,0);
            AV = A*V;

            v = [];
            W = AV;
            if size(block.V,1) > 0
                v = block.AV' * W;
            end
            y = AV' * W;
            block.N = [ block.N, v ; v', y ];
        
            v = [];
            W = block.g5*AV;
            if size(block.V,1) > 0
                v = block.V' * W;
            end
            y = V' * W;
            block.M = [ block.M, v ; v', y ];
        end
        
        block.V = [block.V,V];
        block.AV = [block.AV,AV];
    end
end

%%
function block = extended_hierarchy_update( A, V, block, settings )
    
    block = extend_subspace( A, V, block, settings );
    
    block = clean_up_subspace( block, settings );
    
    block = get_P_Ac( A, block, settings );
    
    if block.level > 1 && settings.coarsegrid
        PV = [];
        if settings.topdown
            PV = block.P'*block.V;
        end
        block.coarser_block = extended_hierarchy_update( block.coarser_A, PV, block.coarser_block, block.coarser_settings );
    end
end


%%
function block = clean_up_subspace( block, settings )
      if(settings.onD)
        [Z,T] = eig(block.N);
      else
        [Z,T] = eig(block.N,block.M);
      end
      block.T = diag(T);

      [~,I] = sort(abs(diag(T)));
      
      Z = Z(:,I(1:min(settings.Nritz,end)));
      [Z,~] = qr(Z,0);
      
      block.V = block.V * Z;
      block.AV = block.AV * Z;
      
      if(~settings.onD)
          block.M = Z' * block.M * Z;
      end
      block.N = Z' * block.N * Z;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive setup algorithm
%
%%
% initial setup
function [ block, settings ] = init_setup( A, settings )
    
    settings.coarsegrid = 1;
    block = get_blocks(A,settings);
    settings.coarsegrid = 0;
    
%     V = randn(size(A,1),settings.Nvec);
%     [V, ~] = qr(V, 0);
    V = zeros(size(A,1),settings.Nvec);
    for i=1:settings.Nvec
        V(i,i)=1;
    end
    
    % apply smoother a few times
    for i = 1:settings.Ninit
        fprintf('level %d -- performed %d initial setup iterations\n',settings.level,i);
        for k = 1:settings.Nvec
            V(:,k) = amg(A,V(:,k),block,settings,0);
        end
    
        [V, ~] = qr(V, 0);
        
        block.V = V;
    end
    
    settings.coarsegrid = 1;
    block = get_P_Ac( A, block, settings );
    
    if block.level > 1
        block = get_coarser_settings( block, settings );
        [block.coarser_block,block.coarser_settings] = init_setup( block.coarser_A, block.coarser_settings );
    end
end


%%
% iterative part
function block = iterative_setup(A, block, settings)

    W = zeros(size(A,1),settings.Nvec);
    for i = 1:settings.Nsetup
        for k = 1:settings.Nvec
            W(:,k) = amg(A,block.V(:,k),block,settings,0);
        end

        [block.V,~] = qr(W,0);
        
        block = hierarchy_update( A, block, settings );

        fprintf('level %d -- performed %d setup iterations\n',settings.level,i);
    
        block.coarser_settings.Nsetup = i;
        if isequal(settings.type, 'simple')
            if block.level > 1
                block.coarser_block = iterative_setup( block.coarser_A, ...
                    block.coarser_block, block.coarser_settings );
            end
        end
    end
end


%%
function block = hierarchy_update( A, block, settings )

    block = get_P_Ac( A, block, settings );
    
    if block.level > 1 && settings.coarsegrid
        block.coarser_block = hierarchy_update( block.coarser_A, ...
            block.coarser_block, block.coarser_settings );
    end
end


%%
function [ block ] = get_P_Ac( A, block, settings )

    k = length(block.coarseindex);
    n = prod(settings.dim);
    N = settings.Nvec;
    
    block.V(1:5,1);
    block.P = spalloc(n,k*N,n*N);
    for j=1:k
        [block.P_dense(block.coarseindex{j},1:N),~] = qr(block.V(block.coarseindex{j},1:N),0);
        [block.P(block.coarseindex{j},N*(j-1)+1:N*j),~] = qr(block.V(block.coarseindex{j},1:N),0);
    end
    block.coarser_A = (block.P)'*(A*(block.P));

end

