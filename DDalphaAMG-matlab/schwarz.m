function [x,r] = schwarz(A,b,r,x,block,settings)

num_color = max(block.color);
p = length(block.color);

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