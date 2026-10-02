clear;

addpath(genpath('/PATH/TO/MOSEK/'));
addpath(genpath('/PATH/TO/YALMIP/'));

params = struct();
params.time_limit = 10800;
params.max_nodes = 1000;
params.opt_gap = 1e-4;
params.sdp_verbose = 0;
params.n_threads = 10;
params.cp_max_iter = 100;
params.cp_eps_ineq = 1e-2;
params.cp_max_sep_ineq = 100000;
params.cp_max_new_tri = 1000;
params.cp_max_new_rlt = 1000;
params.cp_max_new_split = 1000;
params.cp_max_new_penta = 1000;
params.cp_max_new_epta = 1000;
params.cp_max_new_enna = 1000;
params.cp_ntimes_penta = 500;
params.cp_ntimes_epta = 1000;
params.cp_ntimes_enna = 2000;
disp(params)

instances = [
    [80 50 1];
    [80 50 2];
    [80 50 3]
];
num_instances = size(instances, 1);

type = "BW1"; % BW1 - Type1; BW2 - Type2; BW3 - Type3

for t = 1:num_instances

    n = instances(t,1);
    d = instances(t,2);
    s = instances(t,3);
    
    fprintf("\n\t Type: %s \tInstance: %d-%d-%d\n\n", type, n, d, s);
            
    data = load("./instances" + type + "/" + num2str(n) + "_" + num2str(d) + "_" + num2str(s) + ".mat");
    Q = data.Q;
    c = data.c;

    [f_best, x_best] = vns_sum_multistart(Q, c, 100);
    result = solve_ternary_mosek_sum(Q, c, f_best, x_best, params);

end
