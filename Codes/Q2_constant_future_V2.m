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

%% one step ahead


n = length(A) - 1;
n1 = length(B_prime) -1;

%% Finding F and G

d =delay*2;
a1 = [A(1:end-delay), zeros(1, d-1)]; 
na1 = length(a1) - 1;
b1 = [1];
D = [1, zeros(1, na1+d-1)]; % *****here we have q^-1 instead of q^1****

[F, G] = solve_diophantin_general(a1, b1, D);

F = remove_first_zeros(F);

assert(length(F)-1==d-1, "Deg. F is not correct")
assert(length(G)-1== na1-1, "Deg. G is not correct")
assert(F(1)==1, "F is not monic")

%% Finding R amd R_bar

LHS = conv(B_prime, F);
LHS = remove_first_zeros(LHS);
% a11 = conv(1, [1, zeros(1, d-delay+1)]);
% b11 = [1]; 
D1 = conv(LHS, [1, zeros(1, d-delay)]);

% [R, R_bar] = solve_diophantin_general(a11, b11, D1)

% R_bar = remove_first_zeros(R_bar); %%%%
% R = remove_first_zeros(R);

R = D1(1:d-delay+1);
R_bar = D1(d-delay+2:end);

assert(length(R)-1==d-delay, "Deg. R is not correct")
assert(length(R_bar)-1 == na1-2, "Deg. R_bar is not correct")

%% run
rng(50);

main_folder_name = "constant-future";
main_folder = 'images/Q2/' + main_folder_name + "/";
sub_name = "disturbance-";
main_title = sub_name;
close_all = true;

% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
y_star = (-1).^ceil(t/200);
y_star(1) = -1;
input_noise_variance = 0.0;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);
y_star = y_star + input_noise;


% noise and its degree
noise_poly = [1].';
len_noise = length(noise_poly);
noise_variance = 0.0; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
% v = [zeros([1,ceil(num_samples/2)]), 10*ones([1,ceil(num_samples/2)])];
v = zeros([1,num_samples]);


len_desA = 6;
len_desB = 2;

% initial conditions for y
skip_instances = d;
total_parameters = len_desA - 1 + len_desB;
y = [];

y(1:skip_instances) = 0;
u(1:skip_instances) = 0;
u_bar(1:skip_instances) = 0;
% other initial conditions:

theta_real = [A(2:end).'; B_prime.'];

% controller parameters
for i = skip_instances:num_samples-d     
    phi_t = [[-y(i-1:-1:max(i-(length(A) - 1), 1)), zeros(1, 1-(i-(length(A) - 1)))], [u(i-delay:-1:max(i-length(B_prime)-delay+1, 1)), zeros(1, 1-(i-length(B_prime)-delay+1))]].';  % added zeros for initial conditions
    noise_t = [noise(i:-1:max(i-(len_noise-1), 1)), zeros(1, 1-(i-(len_noise-1)))] * noise_poly;
    y(i) = phi_t.' * theta_real + noise_t + B * [v(i:-1:max(i-(length(B)-1), 1)), zeros(1, 1-(i-(length(B) - 1)))].';
    u(i) = (y_star(i+d) - filter_signal(G, y, i) - filter_signal(R_bar(1:end), u, i-1))/(sum(R));
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

