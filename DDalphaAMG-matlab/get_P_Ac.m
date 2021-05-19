function [ block ] = get_P_Ac( A, block, settings )

    k = length(block.coarseindex);
    n = prod(settings.dim);
    N = settings.Nvec;
    
    block.P = spalloc(n,k*N,n*N);
    for j=1:k
        [block.P_dense(block.coarseindex{j},1:N),~] = qr(block.V(block.coarseindex{j},1:N),0);
        [block.P(block.coarseindex{j},N*(j-1)+1:N*j),~] = qr(block.V(block.coarseindex{j},1:N),0);
    end
    block.coarser_A = (block.P)'*(A*(block.P));

end