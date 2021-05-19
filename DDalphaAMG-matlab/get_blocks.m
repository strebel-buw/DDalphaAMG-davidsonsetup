function [block] = get_blocks(A,settings)
% ----------------------------------------------------------------------
% function block = get_blocks(dim,div)
%
% input:
%
% A
% square matrix of size n = prod(dim)
%
% settings
% struct that contains all required parameters to compute the domain
% decomposition
%
%
% output:
%
% block
% struct that contains
%
% block.A
% cell array of subdomain operators
%
% block.color
% aray of colors of each subdomain (for sequential multiplicative
% schwarz)
%
% block.index
% cell array that contains the indices of vertices in each subdomain
%
% block.restrictive_index
% if settings.restrictive == 1, then get_blocks returns a cell array
% that contains the indices of vertices, which don't belong to the overlap
% ----------------------------------------------------------------------

% compute number of vertices in whole domain
dim = settings.dim;
div = settings.div;
griddiv = settings.griddiv;
overlap = settings.overlap;
n = prod(dim);
d = length(dim);

% initialize indices for vertices
index = 1:n;
index = reshape(index,dim);

% compute subdomain indices of 'index'
% initialize subdomain counter
for j = 1:d
    count{j}{1} = j;
    count{j}{2} = 1;
    for i = 1:div(j)
        subbegin{j,i} = floor((i-1)*(dim(j)/div(j)))+1 : floor(i*(dim(j)/div(j)));
    end
end

% compute number of colors
if overlap == 0
    if max(div)>1
        num_color = 2;
    else
        num_color = 1;
    end
else
    num_color = 1;
    color_pow(1) = 1;
    for i = 1:d
        if div(i)>1
            num_color = num_color*2;
        end
        color_pow(i+1) = num_color;
    end
end

% compute subdomain indices
flag = 0;
k = 0;
while flag == 0;
    k = k+1;
    
    % specify input arguments for 'index'
    for j = 1:d
        arguments{j} = subbegin{count{j}{:}};
        current_index(j) = count{j}{2};
    end
    % indices of k-th subdomain
    block.index{k} = index(arguments{:});
    block.index{k} = reshape(block.index{k},numel(block.index{k}),1);
    

    % increase counter
    % if counter = (1,...,1), then flag = 1 -> stop
    flag = 1;
    current = 1;
    while flag == 1 && current <= d
        if count{current}{2} < div(current)
            count{current}{2} = count{current}{2} + 1;
            flag = 0;
        else
            count{current}{2} = 1;
            current = current+1;
        end
    end
    
    % compute color
    if overlap == 0
        block.color(k) = mod(sum(current_index),num_color) + 1;
    else
        indicator = 0;
        for i = 1:d
            indicator = indicator + (current_index(i)-1)*color_pow(i);
            %(2^(i-1));
        end
        block.color(k) = mod(indicator,num_color) + 1;
    end
end

% store indices of subdomains without overlap
if settings.restrictive == 1
    block.restrictive_index = block.index;
end

% compute overlap indices
for j = 1:k
    for i = 1:overlap
        [tmp,block.index{j}] = find(A(block.index{j},:));
    end
    block.index{j} = unique(block.index{j});
end

block.level = settings.level;

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% compute next coarser grid

if settings.coarsegrid == 1
    % compute coarser grid indices of 'index'
    % initialize coarser grid vertice counter
    for j = 1:d
        coarsecount{j}{1} = j;
        coarsecount{j}{2} = 1;
        for i = 1:griddiv(j)
            coarsebegin{j,i}= floor((i-1)*(dim(j)/griddiv(j)))+1 : floor(i*(dim(j)/griddiv(j)));
        end
    end

% compute coarser grid indices
    flag = 0;
    k = 0;
    while flag == 0
        k = k+1;
    
        % specify input arguments for 'index'
        for j = 1:d
            coarsearguments{j} = coarsebegin{coarsecount{j}{:}};
        end
        % indices of k-th projector column
        coarsergridindex{k} = index(coarsearguments{:});
        coarsergridindex{k} = reshape(coarsergridindex{k},numel(coarsergridindex{k}),1);
    

        % increase counter
        % if counter = (1,...,1), then flag = 1 -> stop
        flag = 1;
        current = 1;
        while flag == 1 && current <= d
            if coarsecount{current}{2} < griddiv(current)
                coarsecount{current}{2} = coarsecount{current}{2} + 1;
                flag = 0;
            else
                coarsecount{current}{2} = 1;
                current = current+1;
            end
        end
    end
    
    block.coarseindex = coarsergridindex;
end
   
end

    
