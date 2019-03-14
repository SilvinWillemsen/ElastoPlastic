function [B, C, N, h, bB, bC] = createStringElastoPlastic(c, kappa, L, s0, s1, k)
   
    h = sqrt((c^2 * k^2 + 4 * s1 * k + sqrt((c^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
    N = floor(L/h); % Number of gridpoints
    h = L/N; % Recalculate gridspacing
    
    N = N - 2;

    Dxxxx = (sparse(3:N, 1:N-2, ones(1, N-2), N, N) + ...
            sparse(2:N, 1:N-1, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, 6 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, -4 * ones(1, N-1), N, N) + ...
            sparse(1:N-2, 3:N, ones(1, N-2), N, N));
    Dxx =   (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
            sparse(1:N, 1:N, -2 * ones(1, N), N, N) + ...
            sparse(1:N-1, 2:N, ones(1, N-1), N, N));
%     Dxxxx = Dxx*Dxx;
    N = N - 2;
    B = (2 * eye(N) + c^2 * k^2 / h^2 * Dxx(2:end-1,2:end-1) - kappa^2 * k^2 * Dxxxx(2:end-1,2:end-1) / h^4 + 2 * s1 * k * Dxx(2:end-1,2:end-1) / h^2) / (1 + s0 * k);
    C = -((1 - s0 * k) * eye(N) + 2 * s1 * k * Dxx(2:end-1,2:end-1) / h^2) / (1 + s0 * k);
    
    kOh = (kappa/h^2)^2;
    cOh = (c/h)^2;
    phi = (2 * s1)/(k * h^2);
    bB = (sparse(3:N, 1:N-2, kOh * ones(1, N-2), N, N) + ...
                    sparse(2:N, 1:N-1, (-4 * kOh - cOh - phi) * ones(1, N-1), N, N) + ...
                    sparse(1:N, 1:N, (6*kOh + 2 * cOh - 2/k^2 + 2 * phi) * ones(1, N), N, N) + ...
                    sparse(1:N-1, 2:N, (-4 * kOh - cOh - phi) * ones(1, N-1), N, N) + ...
                    sparse(1:N-2, 3:N, kOh * ones(1, N-2), N, N));
    bC = sparse(2:N, 1:N-1, phi * ones(1, N-1), N, N) + ...
                    sparse(1:N, 1:N, (2/k^2 - 2 * phi) * ones(1, N), N, N) + ...
                    sparse(1:N-1, 2:N, phi * ones(1, N-1), N, N);

    
    N = N + 4;
    
    Dxx2 = zeros(N);
    Dxx2(2:end-1, 2:end-1) = Dxx;
    Dxx = Dxx2;
    
    Dxxxx2 = zeros(N);
    Dxxxx2(2:end-1, 2:end-1) = Dxxxx;
    Dxxxx = Dxxxx2;
    
    B2 = zeros(N);
    B2(3:end-2, 3:end-2) = B;
    B = B2;

    C2 = zeros(N);
    C2(3:end-2, 3:end-2) = C;
    C = C2;
    
    bB2 = zeros(N);
    bB2(3:end-2, 3:end-2) = bB;
    bB = bB2;

    bC2 = zeros(N);
    bC2(3:end-2, 3:end-2) = bC;
    bC = bC2;
end