function [ block ] = get_coarser_settings( block, settings )

block.coarser_settings = settings;
block.coarser_settings.level = settings.level-1;
block.coarser_settings.dim = settings.coarse_dim;
block.coarser_settings.div = settings.coarse_div;
block.coarser_settings.griddiv = settings.coarse_griddiv;
block.coarser_settings.coarse_dim = block.coarser_settings.coarse_griddiv;
block.coarser_settings.coarse_dim(1) = block.coarser_settings.coarse_dim(1)*settings.Nvec;
block.coarser_settings.coarsening = settings.coarse_coarsening;

end

