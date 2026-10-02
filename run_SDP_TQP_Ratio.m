clear;

addpath(genpath('/PATH/TO/MOSEK/'));
addpath(genpath('/PATH/TO/YALMIP/'));


%%%%%%%%%%%%%%%%%% TQP-RATIO SDP APPROACH %%%%%%%%%%%%%%%%%%%%%%

params = struct();
params.time_limit = 10800;
params.max_nodes = 10000;
params.opt_gap = 1e-4;
params.sdp_verbose = 0;
params.n_threads = 10;
params.cp_max_iter = 100;
params.cp_eps_ineq = 1e-6;
params.cp_max_sep_ineq = 100000;
params.cp_max_new_tri = 5000;
params.cp_max_new_rlt = 5000;
params.cp_max_new_split = 5000;
params.cp_max_new_penta = 5000;
params.cp_max_new_epta = 5000;
params.cp_max_new_enna = 5000;
params.cp_ntimes_penta = 500;
params.cp_ntimes_epta = 1000;
params.cp_ntimes_enna = 2000;
disp(params)

instances = [
    [70 25 1];
    [70 25 2];
    [70 25 3]
];
num_instances = size(instances, 1);


sense = 0; % 0 - minimize; 1 maximize;
type = "RATIO";

for t = 1:num_instances

    n = instances(t,1);
    d = instances(t,2);
    s = instances(t,3);
    
    fprintf("\n\t Type: %s \tInstance: %d-%d-%d\n\n", type, n, d, s);

    filename = sprintf('./instancesRATIO/%d_%d_%d.mat', n, d, s);
    data = load(filename).data;

    [f_best, x_best] = vns_ratio_multistart(data.A, data.a, data.a0, data.B, data.b, data.b0, 100, sense);
    result = solve_ternary_mosek_ratio(data, f_best, x_best, params);

end



%%%%%%%%%%%%%%%%%% DINKELBACH (PARAMETRIC QUTO) %%%%%%%%%%%%%%%%%


params = struct();
params.time_limit = 3600;
params.max_nodes = 10000;
params.opt_gap = 1e-4;
params.sdp_verbose = 0;
params.n_threads = 10;
params.cp_max_iter = 100;
params.cp_eps_ineq = 1e-2;
params.cp_max_sep_ineq = 100000;
params.cp_max_new_tri = 5000;
params.cp_max_new_rlt = 5000;
params.cp_max_new_split = 5000;
params.cp_max_new_penta = 5000;
params.cp_max_new_epta = 5000;
params.cp_max_new_enna = 5000;
params.cp_ntimes_penta = 500;
params.cp_ntimes_epta = 1000;
params.cp_ntimes_enna = 2000;
disp(params)

instances = [
    [70 25 1];
    [70 25 2];
    [70 25 3]
];
num_instances = size(instances, 1);

tol = 1e-4; % optimality tolerance
sense = 0; % 0 - minimize; 1 maximize;
type = "RATIO";


for t = 1:num_instances

    n = instances(t,1);
    d = instances(t,2);
    s = instances(t,3);
    
    fprintf("\n\t Type: %s \tInstance: %d-%d-%d\n\n", type, n, d, s);

    filename = sprintf('./instancesRATIO/%d_%d_%d.mat', n, d, s);
    data = load(filename).data;

    [f_best, x_best] = vns_ratio_multistart(data.A, data.a, data.a0, data.B, data.b, data.b0, 100, sense);
    lambda_init = f_best;
    result = compute_dinkelbach(lambda_init, data, tol, params, sense); % last parameter: 0 minimize - 1 maximize

end



