function simple_msd

close all
clear all
clc

% You can change these parameters
N = 20; % Number of masses
M = 0.01; % The mass
ks = 1000; % The spring coefficient
bs = 0; % THe damping coefficient
r = 0.001; % Natural length of the spring
L = 0.2; % Total length of the string


% Initialize position of the mass network
for i = 0:N-1
    x(i+1, :) = [0 -L/(N-1)*i];
end
v = zeros(N,2); % velocity
a = zeros(N,2); % acceleration

% Build stiffnes and damping matrix
KS = zeros(N,N);
BS = zeros(N,N);
for i = 1:N-1
    KS(i, i+1) = ks;
    KS(i+1, i) = ks;  
    BS(i, i+1) = bs;
    BS(i+1, i) = bs;  
end

% Time vector
ts = 1e-3;
tsim = 5;
timespan = 0:ts:tsim;

% Peridoc motion input to the the spring
xin = 0.002*sin(2*pi*7.7162*4*timespan);
xindot = gradient(xin, 1);
xinddot = gradient(xin, 2);

xk = x;
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
        fprintf('Initial tension of the string in x: %dN and in y: %dN\n', ...
                abs(FS(1,1)), abs(FS(1,2)));
    end
    
    % Particle-1 is fixed, don't modify its states!
    a(2:end-1, :) = (FS(2:end-1, :)+FD(2:end-1,:))./M;
    v(2:end-1, :) = v(2:end-1,:) + a(2:end-1,:) .* ts;
    x(2:end-1, :) = x(2:end-1,:) + v(2:end-1,:) .* ts;

    xk = [xk;x];
end

figure
hold on
h = plot(0,0);
ylim([-L 0]);
xlim([-0.02 0.02]);

% plot positions of all particles
for j = 1:length(timespan)
    set(h, 'XData', xk((j-1)*N+1:(j-1)*N+N,1), ...
        'YData', xk((j-1)*N+1:(j-1)*N+N,2))
    pause(0.1);
    drawnow;
end

end

% -------------------------------------------------------------------------
function output = fs(xi,xj,r,k)
if k > 0 % avoid division by zero
    xij = xi-xj;
    xij_hat = xij/norm(xij);
    output = -k*(norm(xij)-r)*xij_hat;
else
    output = [0 0];
end
end

% -------------------------------------------------------------------------
function output = fd(vi,vj,b)
if b > 0 % avoid division by zero
    vij = vi-vj;
    output = -b*vij;
else
    output = 0;
end
end
