function [x,block] = amg(A,b,block,settings,setup_flag)

num_color = max(block.color);
p = length(block.color);
x = zeros(size(b));
r = b;

% ----------------------------------------------------------------------
% schwarz iteration

% pre smoothing
if settings.presmoothing == 1
    if settings.smoother == 0
        [x,~] = gmres(A, r, settings.blockiter, settings.blocktol, settings.smoother_iter);
        r = b - A*x;
    else
%         [x,r] = schwarz(A,b,r,x,block,settings);
        for smoothing_step = 1:settings.smoother_iter
            for current_color = 1:num_color
                v = (zeros(size(b)));
                for i = 1:p
                    if block.color(i) == current_color
                        if settings.restrictive == 1
                            w = (zeros(size(b)));
                            if settings.blocksolver == 1
                                w(block.index{i}) = A(block.index{i},block.index{i})\r(block.index{i});
                            elseif settings.blocksolver == 2
                                w(block.index{i}) = fgmres(A(block.index{i},block.index{i}),r(block.index{i}),...
                                    settings.blockiter,settings.blocktol,1);
                            elseif settings.blocksolver == 3
                                [w(block.index{i}),~] = pcg(A(block.index{i},block.index{i})'*A(block.index{i},block.index{i}),...
                                    A(block.index{i},block.index{i})'*r(block.index{i}), settings.blocktol,settings.blockiter);
                            elseif settings.blocksolver == 4
                                if mod(i,2)== 1
                                    [w(block.index{i}),~] = pcg(-A(block.index{i},block.index{i}),-r(block.index{i}), settings.blocktol,settings.blockiter);
                                else
                                    [w(block.index{i}),~] = pcg(A(block.index{i},block.index{i}),r(block.index{i}), settings.blocktol,settings.blockiter);
                                end
                            end
                            v(block.restrictive_index{i}) = w(block.restrictive_index{i});
                        else
                            if settings.blocksolver == 1
                                v(block.index{i}) = A(block.index{i},block.index{i})\r(block.index{i});
                            elseif settings.blocksolver == 2
                                v(block.index{i}) = fgmres(A(block.index{i},block.index{i}),r(block.index{i}),...
                                    settings.blockiter,settings.blocktol,1);
                            elseif settings.blocksolver == 3
                                [v(block.index{i}),~] = pcg(A(block.index{i},block.index{i})'*A(block.index{i},block.index{i}),...
                                    A(block.index{i},block.index{i})'*r(block.index{i}) ,settings.blocktol,settings.blockiter);
                            elseif settings.blocksolver == 4
                                if mod(i,2)== 1
                                    [v(block.index{i}),~] = pcg(-A(block.index{i},block.index{i}),-r(block.index{i}) ,settings.blocktol,settings.blockiter);
                                else
                                    [v(block.index{i}),~] = pcg(A(block.index{i},block.index{i}),r(block.index{i}) ,settings.blocktol,settings.blockiter);
                                end
                            end
                        end
                    end
                end
                x = x+v;
                if settings.multiplicative == 1
                    r = b - A*x;
                end
            end
        end
    end
end

% coarse grid interpolation
if settings.coarsegrid == 1
    r_sub = (block.P)'*r;
    if block.level == 1
        if settings.coarsesolver == 1
            x = x + block.P*(full(block.coarser_A)\r_sub);
        elseif settings.coarsesolver == 2
            [v,~] = gmres(block.coarser_A,r_sub, ...
                min(settings.coarsemaxiter/settings.coarserestart,size(block.coarser_A,1)), ...
                settings.coarsetol,min(settings.coarsemaxiter,size(block.coarser_A,1)) );
            x = x+block.P*v;
        end
    else
        [v,block.coarser_block] = kcycle(block.coarser_A,r_sub,block.coarser_block,block.coarser_settings,setup_flag);
        x = x + block.P*v;
    end
    
    r = b - A*x;
end

% post smoothing
if settings.postsmoothing == 1
    if settings.smoother == 0
        [v,~] = gmres(A, r, settings.blockiter, settings.blocktol, settings.smoother_iter);
        x = x + v;
    else
%         [x,r] = schwarz(A,b,r,x,block,settings);
        for smoothing_step = 1:settings.smoother_iter
            for current_color = 1:num_color
                v = (zeros(size(b)));
                for i = 1:p
                    if block.color(i) == current_color
                        if settings.restrictive == 1
                            w = (zeros(size(b)));
                            if settings.blocksolver == 1
                                w(block.index{i}) = A(block.index{i},block.index{i})\r(block.index{i});
                            elseif settings.blocksolver == 2
                                w(block.index{i}) = fgmres(A(block.index{i},block.index{i}),r(block.index{i}),...
                                    settings.blockiter,settings.blocktol,1);
                            elseif settings.blocksolver == 3
                                [w(block.index{i}),~] = pcg(A(block.index{i},block.index{i})'*A(block.index{i},block.index{i}),...
                                    A(block.index{i},block.index{i})'*r(block.index{i}), settings.blocktol,settings.blockiter);
                            elseif settings.blocksolver == 4
                                if mod(i,2)== 1
                                    [w(block.index{i}),~] = pcg(-A(block.index{i},block.index{i}),-r(block.index{i}), settings.blocktol,settings.blockiter);
                                else
                                    [w(block.index{i}),~] = pcg(A(block.index{i},block.index{i}),r(block.index{i}), settings.blocktol,settings.blockiter);
                                end
                            end
                            v(block.restrictive_index{i}) = w(block.restrictive_index{i});
                        else
                            if settings.blocksolver == 1
                                v(block.index{i}) = A(block.index{i},block.index{i})\r(block.index{i});
                            elseif settings.blocksolver == 2
                                v(block.index{i}) = fgmres(A(block.index{i},block.index{i}),r(block.index{i}),...
                                    settings.blockiter,settings.blocktol,1);
                            elseif settings.blocksolver == 3
                                [v(block.index{i}),~] = pcg(A(block.index{i},block.index{i})'*A(block.index{i},block.index{i}),...
                                    A(block.index{i},block.index{i})'*r(block.index{i}) ,settings.blocktol,settings.blockiter);
                            elseif settings.blocksolver == 4
                                if mod(i,2)== 1
                                    [v(block.index{i}),~] = pcg(-A(block.index{i},block.index{i}),-r(block.index{i}) ,settings.blocktol,settings.blockiter);
                                else
                                    [v(block.index{i}),~] = pcg(A(block.index{i},block.index{i}),r(block.index{i}) ,settings.blocktol,settings.blockiter);
                                end
                            end
                        end
                    end
                end
                x = x+v;
                if settings.multiplicative == 1
                    r = b - A*x;
                end
            end
        end
    end
end

if setup_flag
    block.V = [block.V, x];
end
