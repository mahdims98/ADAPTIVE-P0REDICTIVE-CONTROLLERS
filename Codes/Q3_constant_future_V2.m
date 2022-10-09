%% discretization of the model
clc; clear;
Ts = 0.75;
delay = 3;

z = tf("z", Ts, 'Variable','z^-1');
sysD = z^-delay * (z+1.5) / ((z-0.25)*(z-0.85));

[numD, denD] = tfdata(sysD, 'v');

% bodeplot(sysC);
% hold on;
% bodeplot(sysD);

B = numD;
A = denD;

% removing zeros of B at the begining
B_prime = remove_first_zeros(B);

d = delay*2;

%% initial conditions
A_estimated = [1,0.1,0.1,0.1, 0.1];
B_estimated = [0, 0, 0, 10, 10]; % first element must be zero

len_desA = length(A_estimated);
len_desB = length(B_estimated);

%% run
rng(50);

main_folder_name = "constant-future";
main_folder = 'images/Q3/' + main_folder_name + "/";
sub_name = "base-";
main_title = sub_name;
close_all = true;

% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
y_star = (-1).^ceil(t/200);
y_star(1) = -1;
input_noise_variance = 0.00;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);
y_star = y_star + input_noise;


% noise and its degree
noise_poly = [1,0,0].';
len_noise = length(noise_poly);
noise_variance = 0.01; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
% v = [zeros([1,ceil(num_samples/2)]), 10*ones([1,ceil(num_samples/2)])];
v = zeros([1,num_samples]);


% initial conditions for y
skip_instances = d;
y = [];

y(1:skip_instances) = 0;
u(1:skip_instances) = 0;
u_bar(1:skip_instances) = 0;
% other initial conditions:

total_parameters = len_desA - 1 + len_desB;
theta_real = [A(2:end).'; B_prime.'];
theta_real_toPlot = [A(2:end), zeros(1, delay), B_prime];
theta_hat_toPlot = zeros([num_samples, total_parameters]);
% RLS initial parameters
rls_solver = RLSClass(100 * eye(total_parameters), 0.01 * ones([total_parameters,1]));

instances_to_estimate_d = 30;
threshold = 0.1; % threshold to determine d
window = 5;

B_prime_estimated = B_estimated(4:end);
d_estimated = length(B_estimated) - length(B_prime_estimated);
% RLS initial parameters
total_parameters = len_desA - 1 + len_desB - d_estimated; 

rls_solver_new = RLSClass(100 * eye(total_parameters), 0.1 * ones([total_parameters,1]));

% controller parameters
beta_thres = 0.01;
u_thres = 2; 


% pr_transfer = tf(R, P, Ts);
% rlocus(z^-(delay) * 1/sysD * pr_transfer)
% 

for i = skip_instances:num_samples-d
    phi_t = [[-y(i-1:-1:max(i-(length(A) - 1), 1)), zeros(1, 1-(i-(length(A) - 1)))],...
        [u(i-delay:-1:max(i-length(B_prime)-delay+1, 1)), zeros(1, 1-(i-length(B_prime)-delay+1))]].'; % for simulations

    phi_t_d = [[-y(i-1:-1:max(i-(len_desA - 1), 1)), zeros(1, 1-(i-(len_desA - 1)))],...
        [u(i-d_estimated:-1:max(i-length(B_prime_estimated)-d_estimated+1, 1)), zeros(1, 1-(i-length(B_prime_estimated)-d_estimated+1))]].';% for estimation
    noise_t = [noise(i:-1:max(i-(len_noise-1), 1)), zeros(1, 1-(i-(len_noise-1)))] * noise_poly;
    y(i) = phi_t.' * theta_real + noise_t + B * [v(i:-1:max(i-(length(B)-1), 1)), zeros(1, 1-(i-(length(B) - 1)))].';

    theta_hat_new = rls_solver_new.update_RLS(y(i), phi_t_d);

    A_estimated = theta_hat_new(1:(len_desA - 1)).';
    
    B_prime_estimated = theta_hat_new(len_desA:total_parameters).';
    theta_hat_toPlot(i, :) = [A_estimated, zeros(1, d_estimated), B_prime_estimated];

        
    % controller design
    d = d_estimated*2;
    a1 = [[1,A_estimated], zeros(1, d-3)]; 
    na1 = length(a1) - 1;
    b1 = [1];
    D = [1, zeros(1, na1+d-1)]; % *****here we have q^-1 instead of q^1****
    
    [F, G] = solve_diophantin_general(a1, b1, D);
    
    F = remove_first_zeros(F);
    
%     assert(length(F)-1==d-3, "Deg. F is not correct")
%     assert(length(G)-1== na1-1, "Deg. G is not correct")


    LHS = conv(B_prime_estimated, F);
    LHS = remove_first_zeros(LHS);
    D1 = conv(LHS, [1, zeros(1, d-d_estimated)]);
    
    R = D1(1:d-d_estimated+1);
    R_bar = D1(d-d_estimated+2:end);
    
%     assert(length(R)-1==d-d_estimated, "Deg. R is not correct")
%     assert(length(R_bar)-1 == na1-2, "Deg. R_bar is not correct")
    

    u(i) = ((y_star(i+d) - filter_signal(G, y, i) - filter_signal(R_bar(1:end), u, i-1))/(sum(R)));

    
    if abs(u(i)) > u_thres
        u(i) = u_thres;
    end
        

end



%% plotters

if ~exist(main_folder, 'dir')
   mkdir(main_folder)
end


figure()
subplot(2,1,1);
plot(t(1:end-d),y, 'DisplayName','Real')
hold on; 
plot(t, y_star-input_noise, "--r", 'DisplayName','Desired')
xlabel("sample number");
title("output");
legend('Location','best');
subplot(2,1,2);
plot(t(1:end-d), u, 'DisplayName','control input')
xlabel("sample number");
title("control input");
saveas(gcf, main_folder + sub_name + "y-u" + '.jpeg')


% Theta
f5 = figure();
f5.Position = [500 50 900 900];
for j = 1:length(theta_hat_toPlot(1,:))
    title_text = "θ_%d";
    subplot(ceil(length(theta_hat_toPlot(1,:))/2),2,j);
    plot(1:num_samples, theta_hat_toPlot(:,j), 'DisplayName','Predicted')
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * theta_real_toPlot(j) , 'DisplayName','Real')
    title(sprintf(title_text, j));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,main_folder + sub_name + "-theta" + '.jpeg')
