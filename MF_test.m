%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   min_{U,V} 1/2 ||UV^{\top} - M||_F^2, rank(M) = r
%                   U \in \mathbb{R}^{nxd} V \in \mathbb{R}^{nxd}
%
%
clear; close all; clc;
n = 5;
r = 2;
kappa_list = [10,10000];
p = 1;
T = 1000;
eta = 0.01;
thresh_up = 1e14; thresh_low = 1e-13;
errors_GD = zeros(length(kappa_list), T);
errors_ScaledGD = zeros(length(kappa_list), T);
val = [];


for i_kappa = 1:length(kappa_list)
    %data generation
    kappa = kappa_list(i_kappa);
    U_seed = sign(rand(n, r) - 0.5);
%     U_seed = randn(n,r);
    [U_star, ~, ~] = svds(U_seed, r);
    V_seed = sign(rand(n, r) - 0.5);
%     V_seed = randn(n,r);
    [V_star, ~, ~] = svds(V_seed, r);
    Omega_seed = rand(n, n);
    
    sigma_star = linspace(1, 1/kappa, r);
    L_star = U_star*diag(sqrt(sigma_star));
    R_star = V_star*diag(sqrt(sigma_star));
    X_star = L_star*R_star';
    Omega = Omega_seed < p;
    Y = Omega.*X_star;

    %% ScaledGD_Random
    L = randn(n,r)/10;
    R = randn(n,r)/10;


    for t = 1:T
        X = L*R';          
        [~,sl,~] = svd(L*R','econ');
        [P,sr,~] = svd(L,'econ');
        val= [val;[(diag(sl))',(diag(sr))']];
        error = norm((X - X_star), 'fro');
        errors_GD(i_kappa, t) = error;
        if ~isfinite(error) || error > thresh_up || error < thresh_low
            break;
        end
        
        % update L
        Z = X - Y;
        L_plus = L - eta/p*Z*R;%/(R'*R));
        Z = L*R'-Y;
        R_plus = R - eta/p*Z'*L;%/(L'*L );

        L = L_plus;
        R = R_plus;
        
        
    end
    

    
    %% AltScaledGD

    L = randn(n,r)/10;
    R = randn(n,r)/10;


    for t = 1:T
        X = L*R';
        error =  norm((X - X_star), 'fro')+norm((X'-X_star'),'fro');
        errors_ScaledGD(i_kappa, t) = error;
        if ~isfinite(error) || error > thresh_up || error < thresh_low
            break;
        end
       

        X = L*R';
        Z = X - Y;
        L = L - eta/p*Z*R/(R'*R));
        X = L*R';
        Z = L*R'-Y;
        R = R - eta/p*Z'*L/(L'*L);

        
    end
 
end




clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,.5,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 20);
lgd = {};
for i_kappa = 1:length(kappa_list)
    kappa = kappa_list(i_kappa);
    errors = errors_GD(i_kappa, :);
    errors = errors(errors > thresh_low);
    t_subs = 1:1:length(errors);
    semilogy(t_subs-1, errors(t_subs), 'Color', clrs{1}, 'Marker', mks{i_kappa}, 'MarkerSize', 9);
    hold on; grid on;
    lgd{end+1} = sprintf('$\\mathrm{GD}~\\kappa=%d$', kappa);
end
for i_kappa = 1:length(kappa_list)
    kappa = kappa_list(i_kappa);
    errors = errors_ScaledGD(i_kappa, :);
    errors = errors(errors > thresh_low);
    t_subs = 1:1:length(errors);
    semilogy(t_subs-1, errors(t_subs), 'Color', clrs{2}, 'Marker', mks{i_kappa}, 'MarkerSize', 9);
    hold on; grid on;
    lgd{end+1} = sprintf('$\\mathrm{AltScaledGD}~\\kappa=%d$', kappa);
end
xlabel('Iteration count');
ylabel('Relative error');
legend(lgd, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 24);
fig_name = sprintf('MC_n=%d_r=%d_p=%g', n, r, p);
