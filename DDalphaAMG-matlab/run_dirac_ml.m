
rng('default');

% ----------------------------------------------------------------------
%% load system
% 
load 4x4x4x4b6.0000id3n1.mat
m0_crit = 0.7867;

% load 8_wilson_dirac_m0_n1
% m0_crit = 0.7972;

% load 16x16x16x16b6.0000id3n1.mat
% m0_crit = 0.809663251530160;
n = size(D,1);
I = speye(size(D));
g5 = kron(speye(n/12),diag([-ones(1,6) ones(1,6)]));
% ----------------------------------------------------------------------
% systemsize
d = sqrt(sqrt(n/12));
dim           = [12,d,d,d,d];
% subdomainsize = [12,2,2,2,2];
subdomainsize = [12,4,4,4,4];

% ----------------------------------------------------------------------
% parameters
kcycle_length = 5;
Nlevels = 2;
Ninit = 1;
Nsetup = 0;
test_vectors = 20;
davidson_MG_update = 20;
setup_tol = 1e-4;
subspace_size = 2*test_vectors;

% ----------------------------------------------------------------------
settings.onD = 0;
settings.topdown = 1;
settings.bottomup = 0;
settings.Nbottomup = 10;
settings.done = 0;

% ----------------------------------------------------------------------
% first coarsening
coarsening = [6,2,2,2,2];
% coarsening = [6,4,4,4,4];

% later coarsenings
% ---------------------------------
coarse_coarsening = [test_vectors,2,2,2,2];
coarse_subdomainsize = [test_vectors*dim(1)/coarsening(1),2,2,2,2];

% ----------------------------------------------------------------------
shift = 0;
A = D - m0_crit*I;
Q = g5*A;
A_shifted = A - shift*I;
Q_shifted = Q - shift*I;
g5Q_shifted = A - shift*g5;

% b = sum(D,2);
b = ones(n,1);

% ----------------------------------------------------------------------
% translate initialized values
d = length(dim);
div = dim./subdomainsize;
griddiv = dim./coarsening;
coarse_dim = griddiv;
coarse_dim(1) = coarse_dim(1)*test_vectors;
coarse_div = coarse_dim./coarse_subdomainsize;
coarse_griddiv = coarse_dim./coarse_coarsening;

% save initialized values
settings.multigrid = 1;
settings.level = Nlevels-1;
settings.Nlevel = Nlevels-1;
settings.Nvec = test_vectors;
settings.Ninit = Ninit;
settings.Nsetup = Nsetup;
settings.Nblocked = davidson_MG_update;
settings.Nritz = subspace_size;
settings.setup_tol = setup_tol;

settings.coarsegrid = 1;
settings.cycle = kcycle_length;
settings.smoother = 0; % 0 = GMRES, 1 = SAP
settings.presmoothing = 0;
settings.postsmoothing = 1;
settings.smoother_iter = 4;
settings.overlap = 0;
settings.restrictive = 0;
settings.multiplicative = 1;

settings.dim = dim;
settings.div = div;
settings.griddiv = griddiv;
settings.coarsening = coarsening;
settings.coarse_dim = coarse_dim;
settings.coarse_div = coarse_div;
settings.coarse_griddiv = coarse_griddiv;
settings.coarse_coarsening = coarse_coarsening;

% ----------------------------------------------------------------------
% solver settings
settings.tol = 1e-8;
settings.maxiter = 700;
settings.restart = 100;
settings.maxrestart = ceil(settings.maxiter/settings.restart);
% ------------------------------------------

settings.coarsesolver = 2; %1 = matlab \ 2 = gmres
settings.coarsemaxiter = 1000;
settings.coarserestart = 1;
settings.coarsetol = 1e-1;

settings.blocksolver = 2; %1 = matlab \ 2 = gmres
settings.blockiter = 8;
settings.blocktol = 1e-8;

% ------------------------------------------
settings.disp = 1;
% figure(77)
if settings.smoother == 0
    titlename = sprintf('%d^4, %d-lvl, %d blockiter GMRES, shift=%f',dim(2), Nlevels, settings.smoother_iter*settings.blockiter,shift);
else 
    titlename = sprintf('%d^4, %d-lvl, %d blockiter SAP, shift=%f',dim(2), Nlevels, settings.blockiter,shift);
end

% ----------------------------------------------------------------------
%% solve Ax=b

setupTypes = {'simple'};
setupBlocking = {'blocked'};

settings.type = setupTypes;
settings.blocked = setupBlocking;
settings.subdomainsize = subdomainsize;

stats = [];
legendText = [];
i=0;
% ----------------------------------------------------------------------


%% evs vs singvals setup D
% range = 0:5;
range = 0;

%Eigenvector setup with Davidson
res2 = zeros(1,length(range));
time2 = zeros(1,length(range));
i = i+1;
stngs{i} = settings;
k=0;
stngs{i}.blocked = 'davidson';
stngs{i}.onD = 0;
stngs{i}.simpleres = 0;
for steps = range
    stngs{i}.Nsetup = steps;
    k = k+1;
    tic
    block = adaptive_setup( A, stngs{i} );
    [~,res] = outer_fgmres(A,b,block,stngs{i},0);
    time2(k) = toc;
    res2(k) = numel(res); hold on;
end
% yyaxis left
% p2 = plot(range,res2, '-g');
% yyaxis right
% plot(range,time2, '--g');
% legendText{i} = 'Davidson on D';
