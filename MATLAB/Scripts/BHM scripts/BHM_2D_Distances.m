L = 3; % 4 6 8 10

N = L^2; coord = zeros(N,2);
for n = 1:N
    coord(n,:) = [1+mod(n-1,L), ceil(n/L)];
end

dx = zeros(N); dy = zeros(N);
for n = 1:N
    for m = 1:N
        dx(n,m) = min([abs(coord(n,1)-coord(m,1)),3-abs(coord(n,1)-coord(m,1))]);
        dy(n,m) = min([abs(coord(n,2)-coord(m,2)),3-abs(coord(n,2)-coord(m,2))]);
    end
end

dist = sqrt(dx.^2 + dy.^2);