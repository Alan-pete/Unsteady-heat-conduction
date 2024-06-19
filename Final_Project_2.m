close all; clear all; clc;

% variables
xmax = 1; % side of square
ymax = 1; % side of square
tmax = 10; % 10 second time limit
N = 5; % determines grid size (resolution)
n = N+1;
alpha = 1;
gamma = 0.25; % CFL number
dx = xmax/N; % change in x (space)
dy = ymax/N; % change in y (space)
dt = (gamma * dx^2) / alpha; % time step

x_list = []; % list for x-values
y_list = []; % list for y-values
t_list = []; % list for time values

x_range = xmax/dx + 1; % range of x-values
y_range = ymax/dy + 1; % range of y-values
t_range = tmax/dt + 1; % range of time values

for i = 1:t_range
    t_list(end+1) = i-1; % adds an integer to the time list
    x_list(end+1) = i-1; % adds an integer to the x-list
end

temp = zeros(n, n, 1); % creates a 3-D matrix for scheme (i,j,k)

% initial conditions
for k = 1:2
    if k > 1 % if k is not the initial time matrix
        temp(:,:,k) = temp(:,:,k-1); % sets the next time matrix equal to the previous time matrix
    end

    for i = 1:n
        for j = 1:n
            temp(i,j,1) = 100; % sets initial temperature values everywhere to 100
            if j == n
                temp(i,j,2) = 200; % right side of square at 200 degrees for t=1
            end
            if i == n
                temp(i,j,2) = 200; % top side of square at 200 degrees for t=1
            end
        end
    end
end

% Unsteady 2-D heat conduction function
conduct = @(T_n_i_j, T_n_ip1_j, T_n_in1_j, T_n_i_jp1, T_n_i_jn1, g) T_n_i_j + gamma*(T_n_ip1_j - 2*T_n_i_j + T_n_in1_j) + gamma*(T_n_i_jp1 - 2*T_n_i_j + T_n_i_jn1);

% numerical scheme (3-nested for loop)
for k = 2:t_range
    if k > 2
        temp(:,:,k) = temp(:,:, k-1); % sets current values to previous values for this time step calculations
    end

    for i = 2:n-1 % keeps bottom row at 100
        for j = 2:n-1 % keeps left column at 100
            % placeholder variables
            a = i-1;
            b = i+1;
            c = j-1;
            d = j+1;
            % function to get new temperature values
            temp_new = conduct(temp(i,j,k-1), temp(b,j,k-1), temp(a,j,k-1), temp(i,d,k-1), temp(i,c,k-1), gamma);
            temp(i,j,k) = temp_new; % puts the calculated value into temperature matrix
            % setting corner values to 150 (median value between 100 and 200)
            temp(1,end,k) = 150;
            temp(end,1,k) = 150;
        end
    end
end

% changing the x, y, and t lists to the correct values for the plots
x_list2 = x_list * dx;
y_list2 = y_list * dy;
t_list2 = t_list * dt;

% Contour plot
figure()
contourf(temp(:,:,end)) % plots the temperature contour at t=10
[contourplot, h] = contourf(temp(:,:,end));
clabel(contourplot, h) % labels the temperature of each contour line
title('Temperature Contour Plot at t=10')
xlabel('X Location / \DeltaT')
ylabel('Y Location / \DeltaT')

% Temperature vs Time plots
figure()
tiledlayout(2,1)

% matrix temperaeture locations
x1 = 0.4/dx;
y1 = 0.4/dy;
x2 = 0.8/dx;
y2 = 0.8/dy;

plot1 = []; % empty list for first plot location
plot2 = []; % empty list for second plot location
for i = 1:t_range
    plot1(i) = temp(x1, y1, i); % puts all the temps at (0.4, 0.4) into the same list
    plot2(i) = temp(x2, y2, i); % puts all the temps at (0.4, 0.4) into the same list
end

% temperature vs time at (0.4,0.4)
nexttile
hold on
plot(t_list2, plot1)
title('Temperature vs Time at (x,y) = (0.4,0.4)')
xlabel('Time (sec)')
ylabel('Temperature (^\circC)')
hold off

% temperature vs time at (0.8,0.8)
nexttile
hold on
plot(t_list2, plot2)
title('Temperature vs Time at (x,y) = (0.8,0.8)')
xlabel('Time (sec)')
ylabel('Temperature (^\circC)')
hold off





