function simple_msd

close all
clear all
clc
N = 20;
M = 0.01;
r = 0.001; % Natural length of the spring
for i = 0:N-1
    x(i+1, :) = [0 -0.2/N*i];
end
v = zeros(N,2); % velocity
a = zeros(N,2); % acceleration

% Stiffness Matrix
ks = 1000;
KS = [0	ks	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0; ...
ks	0	ks	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	ks	0	ks	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	ks	0	ks	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	ks	0	ks	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	ks	0	ks	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	ks	0	ks	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	ks	0	ks	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	ks	0	ks	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	ks	0	ks	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	ks	0	ks	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	ks	0	ks	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	ks	0	ks	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	ks	0	ks	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	ks	0	ks	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	ks	0	ks	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	ks	0	ks	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	ks	0	ks	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	ks	0	ks; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	ks	0];


% Damping Matrix
bs = 0;
BS = [0	bs	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
bs	0	bs	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	bs	0	bs	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	bs	0	bs	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	bs	0	bs	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	bs	0	bs	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	bs	0	bs	0	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	bs	0	bs	0	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	bs	0	bs	0	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	bs	0	bs	0	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	bs	0	bs	0	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	bs	0	bs	0	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	bs	0	bs	0	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	bs	0	bs	0	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	bs	0	bs	0	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	bs	0	bs	0	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	bs	0	bs	0	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	bs	0	bs	0; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	bs	0	bs; ...
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	bs	0];

% Time vector
ts = 1e-3;
tsim = 5;
timespan = 0:ts:tsim;
xk = x;


xin = 0.002*sin(2*pi*7.6948*4*timespan);
xindot = gradient(xin, 1);
xinddot = gradient(xin, 2);

% Perform the simulation
for k = 1:length(timespan)
    x(1,1) = xin(k);
    v(1,1) = xindot(k);
    a(1,1) = xinddot(k);
    FS = zeros(N,2);
    FD = zeros(N,2);
    for i = 1:N
        for j = 1:N
            FS(i,:) = FS(i,:) + fs(x(i,:), x(j,:), r, KS(i,j));
            FD(i,:) = FD(i,:) + fd(v(i,:), v(j,:), BS(i,j));
        end
    end
    
    if k == 1
        disp(FS);
    end
    
    % Particle-1 is fixed, don't modify its states!
    a(2:end-1, :) = (FS(2:end-1, :)+FD(2:end-1,:))./M;
    v(2:end-1, :) = v(2:end-1,:) + a(2:end-1,:) .* ts;
    x(2:end-1, :) = x(2:end-1,:) + v(2:end-1,:) .* ts;
    
    %a(end, 1) = (FS(end, 1)+FD(end,1))./M;
    %v(end, 1) = v(end,1) + a(end,1) .* ts;
    %x(end, 1) = x(end,1) + v(end,1) .* ts;

    xk = [xk;x];
end

figure
hold on
h = plot(0,0);
ylim([-0.2 0]);
xlim([-0.02 0.02]);

% plot positions of all particles
for j = 1:length(timespan)
    set(h, 'XData', xk((j-1)*N+1:(j-1)*N+N,1), 'YData', xk((j-1)*N+1:(j-1)*N+N,2))
    pause(0.1);
    drawnow;
end

end

function output = fs(xi,xj,r,k)
if k > 0 % avoid division by zero
    xij = xi-xj;
    xij_hat = xij/norm(xij);
    output = -k*(norm(xij)-r)*xij_hat;
else
    output = [0 0];
end
end

function output = fd(vi,vj,b)
if b > 0 % avoid division by zero
    vij = vi-vj;
    output = -b*vij;
else
    output = 0;
end
end
