close all; clear all; clc;

% variables
xmax = 10.125; % length of the beam
tmax = 10; % 10 second time limit
alpha = 1;
gamma = 0.5; % CFL number for part A
gamma2 = 1/6; % CFL numberfor part B
dx = 0.25; % change in x (space)
dt = (gamma * dx^2) / alpha; % time step for part A
dt2 = (gamma2 * dx^2) / alpha; % time step for part b
xlim = xmax + dx/2; % makes x_range more divisible by dx (gives whole number)

x_list = []; % list for x values
t_list = []; % list for time values part A
t_list2 = []; % list for time values part B

x_range = xlim/dx + 1;
t_range = tmax/dt + 1; % time range for part A
t_range2 = tmax/dt2 + 1; % time range for part B
temp = zeros(t_range, round(x_range)-1); % matrix for temperature values (part A)
temp2 = zeros(t_range2, round(x_range)-1); % matrix for temperature values (part B)

for i = 1:t_range
    t_list(end+1) = i-1; % adds an integer to the time list part A
end

for i = 1:t_range2
    t_list2(end+1) = i-1; % adds an integer to the time list part B
end

% initial conditions
for i = 1:x_range
    x_list(end+1) = i-1; % adds an integer to the space list
    temp(1,i) = 100; % initial temperature conditions at left end of beam (part A)
    temp2(1,i) = 100; % initial temperature conditions at left end of beam (part B)
end

% Unsteady 1-D heat conduction function
conduct = @(T_n_j, T_n_jp1, T_n_jn1, g) T_n_j + g*(T_n_jp1 - 2*T_n_j + T_n_jn1); % part A
conduct2 = @(T_n_j, T_n_jp1, T_n_jn1, g) T_n_j + g*(T_n_jp1 - 2*T_n_j + T_n_jn1); % part B

% scheme for part A
for i = 2:t_range
    for j = 1:x_range
        a = j-1; % value at far left of the beam
        if a < 1 % outside boundaries
            a = 1;
        end
        b = j+1; % value at far right of the beam
        if b > x_range % outside boundaries
            b = x_range;
        end
        
        % numerical scheme
        temp_new = conduct(temp(i-1,j), temp(i-1,b), temp(i-1,a), gamma);
        temp(i,j) = temp_new;
        if j == x_range
            temp(i,j) = 200; % right end of the beam stays at 200
        end
    end
end

% scheme for part B
for i = 2:t_range2
    for j = 1:x_range
        a = j-1; % value at far left of the beam
        if a < 1 % outside boundaries
            a = 1;
        end
        b = j+1; % value at far right of the beam
        if b > x_range % outside boundaries
            b = x_range;
        end
        
        % numerical scheme
        temp_new2 = conduct(temp2(i-1,j), temp2(i-1,b), temp2(i-1,a), gamma2);
        temp2(i,j) = temp_new2;
        if j == x_range
            temp2(i,j) = 200; % right end of the beam stays at 200
        end
    end
end

% plotting
x_list2 = x_list * dx; % gets the appropriate interval sizes for plot
t_list2A = t_list * dt; % gets the appropriate time step sizes for plot
t_list2B = t_list2 * dt2; % gets the appropriate time step sizes for plot

tiledlayout(3,1)

% plotting temperature versus time at x = 0
nexttile
hold on
plot(t_list2A, temp(:,1), t_list2B, temp2(:,1), 'linewidth', 1)
xlabel('Time (s)')
ylabel('Temperature (^\circC)')
legend('\gamma = 0.5', '\gamma = 1/6', 'location', 'northwest')
title('Temperature vs Time at x = 0')
hold off

% plotting temperature versus time at x = 5.125
nexttile
hold on
plot(t_list2A, temp(:,5.25/dx + 1), t_list2B, temp2(:,5.25/dx + 1), 'linewidth', 1)
xlabel('Time (s)')
ylabel('Temperature (^\circC)')
legend('\gamma = 0.5', '\gamma = 1/6', 'location', 'northwest')
title('Temperature vs Time at x = 5.125')
hold off

% plotting temperature distribution at time t = 10 seconds
nexttile
hold on
plot(x_list2, temp(end,:), x_list2, temp2(end,:), 'linewidth', 1)
xlabel('Position')
ylabel('Temperature (^\circC)')
legend('\gamma = 0.5', '\gamma = 1/6', 'location', 'northwest')
title('Temperature Along the Beam at t = 10 sec')
hold off

