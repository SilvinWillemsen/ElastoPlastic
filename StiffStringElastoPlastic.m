% clear all;
% close all;
clc;

Fs = 44100;     % Sampling rate
f0 = 196.00;    % G3
c = f0 * 2;     % Wave-Speed
gamma = f0*2;     % Wave-Speed
L = 1;
k = 1/Fs;       % Time-step
s0 = 0.1;       % Damping coefficients
s1 = 0.005;

E = 2e11;
r = 0.0005; 
A = r^2*pi;
I = pi*r^4 / 4;
rho = 7850;
kappa = sqrt(E*I/(rho*A));
kappa = 2;
% B = 0.0001; %inharmonicity coefficient
% kappa = sqrt(B)*(gamma/pi); % Stiffness Factor     
% 
% % Calculate grid spacing
% h = sqrt((gamma^2*k^2 + 4 * s1 * k + sqrt((gamma^2 * k^2 + 4 * s1 * k)^2 + 16 * kappa^2 * k^2)) / 2);
% N = floor(1/h); % Number of gridpoints
% h = 1/N; % Recalculate gridspacing
% 
% % Courant numbers
% lambdaSq = (gamma*k/h)^2; 
% muSq = (k * kappa / h^2)^2;

[B, C, N, h, bB, bC] = createStringElastoPlastic(c, kappa, L, s0, s1, k);

% Initialise state vectors
uPrev = zeros(N,1);
u = zeros(N,1);
uNext = zeros(N,1);  

% Raised cosine
% width = 10;
% loc = 0.5;
% startIdx = floor(floor(loc * N) - width / 2);
% endIdx = floor(floor(loc * N) + width / 2);
% 
% u(startIdx : endIdx) = (1 - cos(2 * pi * [0:width] / width)) / 2;
% uPrev = u;

lengthSound = Fs * 5; % Set the length of the output
out = zeros(lengthSound, 1);

% Boundary condition (clamped, ss, free)
bc = "clamped";

%Bow Model
a = 100;                % free parameter
BM = sqrt(2*a)*exp(1/2);

% User variables
Vb = -0.2;               % Bowing speed
Fb = 50;                % Bowing force / total mass of bow
pickup = floor(N/3);    % Pickup position

% Initialise variables for Newton Raphson 
tol = 1e-4;
qSave = zeros (lengthSound, 1);
qPrev = 0;

drawString = true;
            
bp = 0.45; 
I = zeros(N,1);
I(floor(bp * N)) = 1;
J = 1/h * I;
%%%% the Contact Force (be with you) %%%%%%%%%
mus=0.8;               % static friction coeff
mud=0.3;              % dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
strv=1e-1;             % "stribeck" velocity

Fn = 10;

fc=mud*Fn;             % coulomb force
fs=mus*Fn;             % stiction force
ssparam=[fc,fs,strv];            
  
sig0=1e4;           % bristle stiffness
sig1=.1*sqrt(sig0);   % bristle damping
sig2=0.4;            % viscous friction term 
% sig0=4000;           % bristle stiffness
% sig1=0;   % bristle damping
% sig2=0.25;            % viscous friction term 
sigma=[sig0,sig1,sig2];

z_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!) 

z1=0;               % initial mean bristle displacement
y1=0;               % initial non-linear function
f_fr1=0;                % initial friction force 
f_tot_r1=0;         % initial total force on resonator
f_tot_b1=0;         % initial total force on bow

VbInit = -0.1;
a_z=[1,1/(2*Fs)];
zPrev = 0;
z = 0;

%%%%%% Simple Bow Model %%%%%%%
a = 100;                % free parameter
BM = sqrt(2*a)*exp(1/2);

% Initialise variables for Newton Raphson 
tol = 1e-7;
qSave = zeros (lengthSound, 1);
qPrev = 0;
zdot = 0;
zDotPrev = 0;
bowModel = "elastoPlastic";
Vrel = VbInit;
VrelPrev = VbInit;
excitation = 0;

K1 = -sig1 / ((sig2 + 2/k + 2*s0) * h);
for t = 1 : lengthSound
    scalar = sin(10 * pi * t / Fs);
% scalar = 1;
    Vb = VbInit * scalar;
    if bowModel == "simple"
        b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
        eps = 1;
        i = 0;
        while eps>tol
            q=qPrev-(Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0*qPrev+b)/...
             (Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0);
            eps = abs(q-qPrev);
            qPrev = q;
            i = i + 1;
            if i > 10000
                disp('Nope');
            end
        end
        excitation = k^2*J*Fb*BM*q*exp(-a*q^2);
    elseif bowModel == "elastoPlastic"
        
        % calculate pre-computable part of the FDS
        b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
        if fc>0
            eps = 1;
            i = 0;
            while eps>tol && i < 100
                espon=exp(-(Vrel/strv).^2);         %exponential function
                zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
                if Vrel==0
                  zss=fs/sig0;
                end
%                 zss = abs(zss);
                %elasto-plastic function \alpha (v,z)
                zssPrev = zss;
                zss = abs(zss);
                alpha = 0;
                if (sign(z)==sign(Vrel))
                    if ((abs(z)>z_ba) && (abs(z)<zss))
                        arg=pi*(z-sign(z)*0.5*(zss+z_ba))/(zss-z_ba);
                        sinarg=sin(sign(z)*arg);
                        alpha=0.5*(1+sinarg);
                    elseif (abs(z)>zss)
                        alpha=1;
                    end
                end
                zss = zssPrev;
                %non-linear function estimate
                fnl=Vrel*(1-alpha*z/zss);
                %% compute derivatives
                
                %dz_ss/dv
                dz_ss = (-2*Vrel*sign(Vrel) / (strv^2*sig0))*(fs-fc)*espon;
                
                dalpha_v=0; %d(alpha)/dv 
                dalpha_z=0; %d(alpha)/dz
                
                zss = abs(zss);
                if ((sign(z)==sign(Vrel)) && (abs(z)>z_ba) && (abs(z)<zss) )
                    cosarg=cos(arg);
                    dalpha_v=0.5*pi*cosarg*dz_ss*(z_ba-z)/(zss-z_ba)^2; 
                    dalpha_z=0.5*pi*cosarg/(zss-z_ba);
                end
                zss = zssPrev;
                
                d_fnlv = 1-z * ((alpha +Vrel*dalpha_v)*zss -dz_ss*alpha*Vrel)/zss^2;
                d_fnlz = -Vrel/zss*(z*dalpha_z +alpha);
                K1 = -k^2 / (h * (1 + s0 * k));
                d_fnl = K1 * d_fnlv + d_fnlz * k / 2;
                
                zdotNext = zdot - (fnl - zdot)/(d_fnl - 1);
                eps = abs(zdotNext-zdot);
                zdot = zdotNext;
                
                z = zPrev + k / 2 * zDotPrev + k / 2 * zdot;
                Vrel = (-sig0 * z - sig1 * zdot - b) / (sig2 + 2/k + 2*s0);
                i = i + 1;
            end
            F = sig0 * z + sig1 * zdot + sig2 * Vrel;
        else 
            zdot=0; err=0; count=0;
            z=0; f_fr=0; v=vs;
        end
        excitation = k^2*J*F;
    elseif bowModel == "hyperbolic"
        v = 1/k * (I' * u - I' * uPrev);
        mu = (mud + (mus - mud) * VbInit/2) / (VbInit/2 + v - Vb);
        excitation = k^2*J*Fn*mu;
    end
    uNext = (B * u + C * uPrev) - excitation / (1 + s0 * k);
    zSave(t) = z;
    zDotSave(t) = zdot;
    zPrev = z;
    zDotPrev = zdot;
    % Plot
    if drawString == true && mod(t,100) == 0 %&& t > Fs * 2
        clf;
        subplot(3,1,1)
        plot(uNext);
        subplot(3,1,2)
        plot(zSave(1:t))
        subplot(3,1,3)
        plot(zDotSave(1:t))
%         ylim([-1e-6 1e-6])
%         hold on;
%         scatter(floor(bp*N), 0);
        
        drawnow;
    end
    alphaSave(t) = alpha;
    % Save output
    out(t) = uNext(pickup);
    
    % update state vectors
    uPrev = u;
    u = uNext;
end

plot(out);