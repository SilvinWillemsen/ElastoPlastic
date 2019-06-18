clear all;
close all;
clc;

%% Drawing functions
drawString = false;
drawStart = 0;
drawspeed = 100;
drawEnergy = false;
hysteresis = false;

%% Choose bow model (elastoPlastic, simple, hyperbolic, cos (raised cosine)
bowModel = "elastoPlastic";
if bowModel ~= "elastoPlastic"
    drawString = false;
end
Fs = 44100;     % Sampling rate
f0 = 196.00;    % G3
c = f0 * 2;     % Wave-Speed
L = 1;          % String length
k = 1/Fs;       % Time-step
s0 = 0.1;       % Damping coefficients
s1 = 0.005;

E = 2e11;       % Young's modulus
r = 0.0005;     % Radius
A = r^2 * pi;   % Crossectional area
Iner = pi*r^4 / 4; % Area moment of inertia
rho = 7850;         % Density
T = c^2 * rho * A;  % Tension
kappa = sqrt(E*Iner/(rho*A)); % Stiffness

%% Create string
[B, C, N, h, Dxx, Dxxxx, s0, s1, bB, bC] = unscaledCreateStringNR(rho, A, T, E, Iner, L, s0, s1, k);

%% Initialise state vectors
uPrev = zeros(N,1);
u = zeros(N,1);
uNext = zeros(N,1);

%% Raised cosine
width = 6;
loc = 0.25;
startIdx = floor(floor(loc * N) - width / 2);
endIdx = floor(floor(loc * N) + width / 2);
if bowModel == "cos"
    u(startIdx : endIdx) = (1 - cos(2 * pi * [0:width] / width)) / 2;
end
uPrev = u;

lengthSound = Fs*2; % Set the length of the output
out = zeros(lengthSound, 1);

%% Initialise energy vectors
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCdamp0StringEnergy = zeros(lengthSound, 1);
rOCdamp1StringEnergy = zeros(lengthSound, 1);
rOCbowStringEnergy = zeros(lengthSound, 1);
rOCenergy = zeros(lengthSound, 1);

pickup = floor(2*N/3);    % Pickup position

bp = 1/4;
bP = floor(bp * N);
I = zeros(N,1);
I(floor(bp * N)) = 1;
J = 1/h * I;

if bowModel == "elastoPlastic"
    %%%% the Contact Force (be with you) %%%%%%%%%
    mus=0.8;               % static friction coeff
    mud=0.3;              % dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
    strv=1e-1;             % "stribeck" velocity
    
    FnInit = 10;
    Fn = FnInit;
    fc=mud*Fn;             % coulomb force
    fs=mus*Fn;             % stiction force
    
    sig0 = 1e4;           % bristle stiffness
    sig1 = 0.1*sqrt(sig0);   % bristle damping
    sig2 = 0.4;            % viscous friction term
    sig3 = 0;
    
    w = (2 * rand(lengthSound,1) - 1);
    
    sigma=[sig0,sig1,sig2];
    
    z_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!)
    
    VbInit = -0.1;
    zPrev = 0;
    z = 0;
    tol = 1e-7;
    
    zdot = VbInit;
    zDotPrev = VbInit;
    an = zdot;
    anPrev = zDotPrev;
    
    Vrel = -VbInit;
    VrelPrev = -VbInit;
    VbPrev = VbInit;
    
    zVec = -z_ba*5:1e-7:z_ba*5;
    vRelVec = -0.5:0.001:0.5;
    for ii = 1 : length(vRelVec)
        espon=exp(-(vRelVec(ii)/strv).^2);%exponential function
        zssVec(ii) =sign(vRelVec(ii)).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
        if vRelVec(ii)==0
            zssVec(ii)=fs/sig0;
        end
    end
    zSave = zeros(lengthSound, 1);
    vRelSave = zeros(lengthSound, 1);
    zDotSave = zeros(lengthSound, 1);
    Fsave = zeros(lengthSound, 1);
    iSave = zeros(lengthSound, 1);
    
else
    a = 100;                % free parameter
    Vb = 0.1;               % Bowing speed
    Fb = 1;                % Bowing force / total mass of string
    BM = sqrt(2*a)*exp(1/2);
    tol = 1e-4;
    
    qSave = zeros (lengthSound, 1);
    qPrev = 0;
end

excitation = 0;

bowVertPosPrev = 0; % for drawing

scalar = 0;


% w = zeros(lengthSound,1);
% wNoise = awgn(w,10);
% plot(wNoise);
zssVecPlotVal = 1;
% figure('Renderer', 'painters', 'Position', [100 100 900 500])
eVec = 2:N-1;
vec = 3:N-2;

for t = 1 : lengthSound
    %     Fn =  t / lengthSound * FnInit
    if bowModel == "elastoPlastic"
        fc=mud*Fn;             % coulomb force
        fs=mus*Fn;             % stiction force
        z_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!)
    end
%     if t < ramp
%         Vb = VbInit * cos(10 * pi * t / Fs);
%     else
        Vb = VbInit;
%     end

    if bowModel == "simple"
        b = 2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev;
        eps = 1;
        i = 0;
        while eps>tol
            q=qPrev-(1/h*Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0*qPrev+b)/...
                (1/h * Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0);
            eps = abs(q-qPrev);
            qPrev = q;
            i = i + 1;
            if i > 10000
                disp('Nope');
            end
        end
        qSave(t) = q;
        Fsave(t) = Fb * sqrt(2 * a) * q * exp(-a*q^2 + 1/2);
        excitation = J * Fsave(t);
        vRelSave(t) = q;
    elseif bowModel == "elastoPlastic"
        
        % calculate pre-computable part of the FDS
        b = (2/k * Vb + 2 * s0 * Vb + I' * bB * u + I' * bC * uPrev);
        if fc>0
            eps = 1;
            i = 0;
            while eps>tol && i < 200
                espon=exp(-(Vrel/strv).^2);         %exponential function
                zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
                if Vrel==0
                    zss=fs/sig0;
                end
                zssSave(t) = zss;
                % elasto-plastic function \alpha (v,z)
                zssPrev = zss;
                zss = abs(zss);
                alpha = 0;

                if (sign(z)==sign(Vrel))
                    if ((abs(z)>z_ba) && (abs(z)<zss))
                        arg=pi*(z-sign(z)*0.5*(zss+z_ba))/(zss-z_ba);
                        sinarg=sin(sign(z)*arg);
                        alpha=0.5*(1+sinarg);
                    elseif (abs(z)>=zss)
                        alpha=1;
                    end
                end
                alphaSave(t) = alpha;
                zss = zssPrev;
                
                fnl = Vrel * (1 - alpha * z / zss);
                
                %% compute derivatives
                
                % Derivative of the steady state function (dz_ss/dv)
                dz_ss = (-2*abs(Vrel) / (strv^2*sig0))*(fs-fc)*espon;
                dz_ssAbs = zss / abs(zss) * dz_ss;
                dalpha_v=0; %d(alpha)/dv
                dalpha_z=0; %d(alpha)/dz
                zss = abs(zss);
                if ((sign(z)==sign(Vrel)) && (abs(z)>z_ba) && (abs(z)<zss))
                    cosarg=cos(sign(z)*arg);
                    dalpha_v=0.5*pi*cosarg*dz_ssAbs*(z_ba-abs(z))/(zss-z_ba)^2;
                    dalpha_z=0.5*pi*cosarg*sign(z)/(zss-z_ba);
                end
                zss = zssPrev;
               
                dfnl_z = -Vrel/zss*(z*dalpha_z + alpha);
                dfnl_v = 1-z * ((alpha+Vrel*dalpha_v)*zss -dz_ss*alpha*Vrel)/zss^2;
                
                K1 = -(2/k + 2*s0 + sig2/(rho * A * h)) / (sig1/(rho * A * h));
                K2 = 2/k;
                dfnl = K1 / 2 * dfnl_v + K2 * dfnl_z - 1;
                
                zdotPrevIt = zdot;
                zdot = zdot - (fnl - zdot) / dfnl;
                z = zPrev + k/2 * (zdot + zDotPrev);
                Vrel = ((-sig0 * z - sig1 * zdot - sig3 * w(t)) / (rho * A * h) - b) / (2/k + 2 * s0 + sig2 / (rho * A * h));
                
                eps = abs(zdot-zdotPrevIt);
                
                i = i + 1;
                
            end
            iSave(t) = i;
            if t > drawStart
                
            end
%             espon=exp(-(Vrel/strv).^2);         %exponential function
%             zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
%             if Vrel==0
%                 zss=fs/sig0;
%             end
%             zssSave(t) = zss;
%             % elasto-plastic function \alpha (v,z)
%             zssPrev = zss;
%             zss = abs(zss);
%             alpha = 0;
% 
%             if (sign(z)==sign(Vrel))
%                 if ((abs(z)>z_ba) && (abs(z)<zss))
%                     arg=pi*(z-sign(z)*0.5*(zss+z_ba))/(zss-z_ba);
%                     sinarg=sin(sign(z)*arg);
%                     alpha=0.5*(1+sinarg);
%                 elseif (abs(z)>=zss)
%                     alpha=1;
%                 end
%             end
%             alphaSave(t) = alpha;
%             zss = zssPrev;
%             zdotPrevIt = zdot;
%             zdot= Vrel * (1 - alpha * z / zss);
%             zdotPrevIt - zdot
            F = (sig0 * z + sig1 * zdot + sig2 * Vrel + sig3 * w(t));
            Fsave(t) = F;
        else
            zdot=0; err=0; count=0;
            z=0; f_fr=0; v=vs;
        end
        excitation = J*F;
        vRelSave(t) = Vrel;
        zSave(t) = z;
        zDotSave(t) = zdot;
        
        % update z and zdot
        zPrev = z;
        zDotPrev = zdot;
        anPrev = an;
        
    elseif bowModel == "hyperbolic"
        v = 1/k * (I' * u - I' * uPrev);
        mu = (mud + (mus - mud) * Vb/2) / (Vb/2 + v - Vb);
        excitation = J*Fn*mu;
    else
        excitation = 0;
    end
    
    %% Update FDS
    uNext = (B * u + C * uPrev) - excitation / (rho * A / k^2 + s0 / k);

    %% Calculate energy of the string
    kinEnergy(t) = rho * A / 2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(t) = T / 2 * 1/h * sum((u(3:N) - u(2:N-1)) .* (uPrev(3:N) - uPrev(2:N-1)))...
        + E * Iner / 2 * 1/h^3 * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1)) ...
        .* (uPrev(eVec+1) - 2 * uPrev(eVec) + uPrev(eVec-1)));
    totEnergy(t) = kinEnergy(t) + potEnergy(t);
    
    
    %% Calculate Rate of Changes of Energy
    rOCkinEnergy(t) = h * rho * A / (2 * k^3) * sum((uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(t) = h * T / (2*k*h^2) * sum((u(vec+1) - 2 * u(vec) + u(vec-1)).* (uNext(vec) - uPrev(vec))) ...
        - h * E * Iner / (2 * k * h^4) * sum((u(vec+2) - 4 * u(vec+1) + 6 * u(vec) - 4 * u(vec-1) + u(vec-2)) .* (uNext(vec) - uPrev(vec)));%...
    rOCdamp0StringEnergy(t) = -2 * s0 * h / (4 * k^2) * sum((uNext - uPrev).*(uNext - uPrev));
    rOCdamp1StringEnergy(t) = 2 * h * s1 / (2 * k^2 * h^2) * sum((u(eVec+1) - 2 * u(eVec) + u(eVec-1) - uPrev(eVec+1) + 2 * uPrev(eVec) - uPrev(eVec-1)) .* (uNext(eVec) - uPrev(eVec)));
    rOCbowStringEnergy(t) = -h * sum(I .* excitation) * (uNext(bP) - uPrev(bP)) / (2 * k);
    rOCenergy(t) = rOCkinEnergy(t) - rOCpotEnergy(t) - rOCdamp0StringEnergy(t) - rOCdamp1StringEnergy(t) - rOCbowStringEnergy(t);
    
    %% Drawing functions
    if drawString == true && mod(t,drawspeed) == 0 && t > drawStart
        if drawEnergy == false 
            clf
            subplot(2,2,1)
            plot(uNext);
           
            %         scatter(repmat(floor(bp*N), 20, 1), [-1e-5:1e-6:1e-5-1e-6]+mod(bowVertPos/20,1e-6),'.');
            text(floor(bp*N) - 2, 0, "$V_B =$ " + num2str(Vb, 2), 'interpreter', 'latex', 'horizontalAlignment', 'right');
            title("String displacement at sample " + num2str(t), 'interpreter', 'latex','Fontsize', 16)
            subplot(2,2,2)
            plot(out(1:t))
            if hysteresis == false
                %% Steady State Curve
                subplot(2,4,5)
                plot(vRelVec, zssVec)
                hold on;
                zssVecVal = find (round(vRelVec * 1e3) == round(Vrel * 1e3));
                if length(zssVecVal) == 1
                    zssVecPlotVal = zssVecVal;
                end
                scatter(round(Vrel * 1e3)*1e-3, zssVec(zssVecPlotVal));
                text(round(Vrel * 1e3)*1e-3, zssVec(zssVecPlotVal) - 1.5e-4, "$v =$ " + num2str(Vrel, 2), 'interpreter', 'latex', 'horizontalAlignment', 'center');
                xlim([-0.5 0.5])
                %         title("$z$", 'interpreter', 'latex')
                xlabel('$v$','interpreter', 'latex')
                ylabel("$z_{ss}(v)$", 'interpreter', 'latex')
                title('Steady-state curve $z_{ss}(v)$', 'interpreter', 'latex', 'Fontsize', 16)
                
                %% Adhesion Map
                subplot(2,4,6)
                alphaPlot = zeros(length(zVec),1);
                zssPlot = abs(zss);
                for ii = 1:length(zVec)
                    if sign(zVec(ii)) == sign(zss)
                        if ((abs(zVec(ii))>z_ba) && (abs(zVec(ii))<zssPlot))
                            arg=pi*(zVec(ii)-sign(zVec(ii))*0.5*(zssPlot+z_ba))/(zssPlot-z_ba);
                            sinarg=sin(sign(zVec(ii))*arg);
                            alphaPlot(ii)=0.5*(1+sinarg);
                        elseif (abs(zVec(ii))>=zssPlot)
                            alphaPlot(ii)=1;
                        end
                    end
                end
                plot(zVec, alphaPlot);
                text(zss + 5e-5, 0.5, '$z_{ss}(v)$', 'interpreter', 'latex', 'horizontalAlignment', 'center')
                
                hold on;
                plot([zss zss], [min(ylim) max(ylim)], '--');
                alphaIdx = find(round(zVec*1e7) == floor(1e7*z));
                if length(alphaIdx) == 1
                    scatter(z, alpha);
                    text(floor(1e5*z)*1e-5+3e-5, alphaPlot(alphaIdx) - 0.05, '$z$', 'interpreter', 'latex', 'horizontalAlignment', 'center');
                end
                %         title('$\alpha$','interpreter', 'latex', 'Fontsize', 18)
                xlabel('$z$','interpreter', 'latex')
                ylabel('$\alpha(v,z)$','interpreter', 'latex')
                title('Adhesion map $\alpha(v,z)$', 'interpreter', 'latex', 'Fontsize', 16)
                
                %% Velocity of Mean Bristle Displacement
                subplot(2,4,7)
                plot(zDotSave(1:t))
                xlabel('$n$ (samples)','interpreter', 'latex')
                ylabel("$\dot z$", 'interpreter', 'latex')
                title('Velocity of the mean bristle displacement', 'interpreter', 'latex', 'Fontsize', 16)
                
                %% Mean Bristle Displacement
                subplot(2,4,8)
                plot(zSave(1:t))
                title("$z$", 'interpreter', 'latex')
                xlabel('$n$ (samples)','interpreter', 'latex')
                ylabel("$z$", 'interpreter', 'latex')
                title('Mean bristle displacement $z$', 'interpreter', 'latex', 'Fontsize', 16)
            else
                subplot(2,1,2)
                numDots = 500;
                mat = [[1:-1/numDots:1/numDots]', [1:-1/numDots:1/numDots]', [1:-1/numDots:1/numDots]'];
                %             vec = t-numDots + 1 : t;
                cla;
                if t > numDots
                                   plot(vRelSave(t-numDots+1:t), Fsave(t-numDots+1:t), 'k', 'Linewidth', 2);
                               hold on;
                    scatter(vRelSave(t-numDots+1:t), Fsave(t-numDots+1:t), 40, mat);
                end
                ylabel("$f(v,z)$", 'interpreter', 'latex')
                xlabel("$v$", 'interpreter', 'latex')
                title('Force function $f(v,z)$ against relative velocity', 'interpreter', 'latex', 'Fontsize', 16)
                %             xlim([-0.25 0.05])
                %             ylim([-0.2 1])
                set(gca, 'Fontsize', 20)
                %             set(gca, 'Position', [0.08 0.11 0.89 0.82])
                grid on;
            end
            drawnow;
            
        elseif drawEnergy == true && t > drawStart - 1000 && mod(t, drawspeed) == 0
            subplot(3,1,1);
            plot(u);
            title("String (state)")
            subplot(3,1,2);
            plot(totEnergy(10:t) / totEnergy(10) - 1);
            title("Total energy")
            subplot(3,1,3);
            cla
            plot(rOCenergy(10:t));
            title("Rate of change of energy + bowing energy (should be 0)");
            %         hold on;
            %         plot(rOCbowStringEnergy(10:t))
            drawnow;
        end
    end
    %     alphaSave(t) = alpha;
    % Save output
    out(t) = uNext(pickup);
    %     out2(t) = u(pickup) + 1 / 2 * (uNext(pickup) - uPrev(pickup));
    
    % update state vectors
    uPrev = u;
    u = uNext;
end

plot(out);