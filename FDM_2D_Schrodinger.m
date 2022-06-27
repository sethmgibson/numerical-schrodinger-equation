% Constants
hbar = 1; % Planck's
m = 1; % Particle mass

% Domain
L = 1;
N = 50;
Nt = 1e3;
dx = 1/N;
dt = 1/Nt;
x = linspace(0, L, N)';
y = linspace(0, L, N);

% Initial Wavefunction
% u0 = sqrt(2).*(sin(pi*x) + sin(pi*y));
u0 = exp(1i*(4*pi*(x+y) + pi/2));

% Potential Function
V = 1e4-1e4.*exp(-(x - L/2).^2/(2*(1/50))).*exp(-(y - L/2).^2/(2*(1/50))); % 2D Gaussian

% V = zeros(N,N); % 2D Square Well
% for j = 1:N/3
%     V(j,:) = 1e7;
%     V(:,j) = 1e7;
%     V(N - j+1,:) = 1e7;
%     V(:,N - j+1) = 1e7;
% end

surf(x,y,V);

% Assembly Matrix
a = dt.*V/(2*hbar);
b = dt*hbar/(4*m*dx*dx);

u0 = reshape(u0, [N*N, 1]);
u = zeros(N*N, Nt);
u(:, 1) = u0;

I = zeros(N, N);
I(1, 1) = 1;
I(N, N) = 1;
A = zeros(N, N);
A(1, 1) = a(1,1) + 4*b;
A(1, 2) = -b;
A(N, N-1) = -b;
A(N, N) = a(N,N) + 4*b;

for j = 2:N-1
    A(j, j-1) = -b;
    A(j, j) = a(j,j) + 2*b;
    A(j, j+1) = -b;

    I(j, j) = 1;
end

% Kronecker Product
B = kron(A, I) + kron(I, A);

% Unitary Matrix
U = (kron(I, I) - 1i.*B)/(kron(I, I) + 1i.*B);

for t = 1:Nt
    u(:, t+1) = U*(u(:, t)./(sqrt(abs(u(:, t)'*u(:, t)))));
end

for t = 1:Nt/6
    u1 = u(:,t);
%     plot(x, real(u(:, t)), x, imag(u(:, t)));
%     plot(x, abs(u(:, t)).^2, 'b', x, V(:)/1e4, 'k', x, real(u(:, t)), 'm', x, imag(u(:, t)), 'c');
    u1 = reshape(u1, [N, N]);
    surf(x, y, abs(u1).^2);
zlim([0 0.05])
    hold on
%     surf(x, y, V/1e6, 'FaceColor','r', 'FaceAlpha',0.5)
    hold off
    pause(0.01)
end