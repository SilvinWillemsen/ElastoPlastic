clear all;
close all;
clc;

Fs = 44100;     % Sampling rate

%% Drawing functions
drawString = true;
drawStart = 0;
drawspeed = 100;
drawEnergy = false;
hysteresis = true;

%% Choose bow model (elastoPlastic, simple, hyperbolic, cos (raised cosine))
bowModel = "elastoPlastic";
if bowModel ~= "elastoPlastic"
%     drawString = false;
end


f0 = 196.00;    % G3
c = f0 * 2;     % Wave-Speed
L = 1;          % String length
k = 1/Fs;       % Time-step
s0 = 0.0;       % Damping coefficients
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
if bowModel == "cos"
    width = 6;
    loc = 0.25;
    startIdx = floor(floor(loc * N) - width / 2);
    endIdx = floor(floor(loc * N) + width / 2);
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
rOCbowEnergyInput = zeros(lengthSound, 1);
rOCbristleEnergy = zeros(lengthSound, 1);

rOCbowEnergy = zeros(lengthSound, 1);
rOCdispEnergy = zeros(lengthSound,1);
rOCbowStringEnergy = zeros(lengthSound, 1);
rOCenergy = zeros(lengthSound, 1);

pickup1 = floor(2*N/3);    % Pickup position
pickup2 = floor(N/2);

bp = 33/N;
bP = floor(bp * N);
I = zeros(N,1);
alph = bp * N - bP;
interPolVec = [(alph * (alph - 1) * (alph - 2)) / -6.0;
    ((alph - 1) * (alph + 1) * (alph - 2)) / 2.0;
    (alph * (alph + 1) * (alph - 2)) / -2.0;
    (alph * (alph + 1) * (alph - 1)) / 6.0];

I(bP-1:bP+2) = interPolVec;
J = 1/h * I;

if bowModel == "elastoPlastic"
    %%%% the Contact Force (be with you) %%%%%%%%%
    mus=0.8;               % static friction coeff
    mud=0.3;              % dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
    strv=1e-1;             % "stribeck" velocity
    
    FnInit = 1;
    Fn = FnInit;
    fc=mud*Fn;             % coulomb force
    fs=mus*Fn;             % stiction force
    
    sig0 = 1e4;           % bristle stiffness
    %     sig1 = 0.001*sqrt(sig0);   % bristle dampin
    sig2 = 0.4;            % viscous friction term
    sig3 = 0.0;
    w = (2 * rand(lengthSound,1) - 1);
    
    z_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!)
    
    VbInit = 0.1;
    zPrev = 0;
    z = 0;
    tol = 1e-7;
    
    zdot = VbInit;
    zDotPrev = VbInit;
    an = zdot;
    anPrev = zDotPrev;
    
%     zdot = 0;
%     zDotPrev = 0;
%     an = 0;
%     anPrev = 0;
%     
    Vrel = -VbInit;
    VrelPrev = -VbInit;
    VbPrev = VbInit;
    %     Vrel = 0;
    %     VrelPrev = 0;
    %     VbPrev = 0;
    
    % for drawing
    zVec = -z_ba*5:z_ba/10:z_ba*5;
    vRelVec = -0.5:0.001:0.5;
    for ii = 1 : length(vRelVec)
        espon=exp(-(vRelVec(ii)/strv).^2);%exponential function
        zssVec(ii) =sign(vRelVec(ii)).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
        if vRelVec(ii)==0
            zssVec(ii)=fs/sig0;
        end
    end
%         sig1 = 4 * sig0 * zssVec(end) / (vRelVec(end) * 1);
    sig1 = 0;
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

scalar = 0;

zssVecPlotVal = 1;
eVec = 2:N-1;
vec = 3:N-2;
saveAlphaItFlag = false;
drawNR = false;
ramp = lengthSound / 2;
figure('Renderer', 'painters', 'Position', [100 100 800 350])

sig1coeff = 3;
for t = 1 : lengthSound
    
    %%% For debugging
    if t == 100
        %         saveAlphaItFlag = true;
    end
    if saveAlphaItFlag
        alphaIt = [];
        vrelIt = [];
        zIt = [];
        zDotIt = [];
        epsSave = [];
    end
    %%%
    %     if t < ramp
    %         Fn = t / ramp * (FnInit - 0.01) + 0.01;
    %     else
    
    %     end
    % FnSave(t) = Fn;
    if bowModel == "elastoPlastic"
        Fn = FnInit;
        fc=mud*Fn;             % coulomb force
        fs=mus*Fn;             % stiction force
        z_ba=0.7*fc/sig0;      % break-away displacement (has to be < f_c/sigma_0!!)
        Vb = VbInit;
    end

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
            VrelPrevIt = tol;
            zPrevIt = tol;
            z_ba=0.7*fc/sig0;
            while eps>tol && i < 50
                espon=exp(-(Vrel/strv).^2);         %exponential function
                zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
                if Vrel==0
                    zss=fs/sig0;
                end
                zssSave(t) = zss;
                zssPrev = zss;
                zss = abs(zss);
                alpha = 0;
                
                % elasto-plastic function \alpha (v,z)
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
                
                an = 2/k * (z - zPrev) - anPrev;
                
                %                 zdotPrevIt = zdot;
                
                zdot = Vrel * (1 - alpha * z / zss);
                
                %                 sig1 = sig1coeff * sig0 * zss * alpha / Vrel;
                
                % functions to perform newton raphson on
                g1 = (2/k + 2 * s0) * Vrel + (sig0 * z + sig1 * zdot ...
                    + sig2 * Vrel + sig3 * w(t)) / (rho * A * h) + b;
                g2 = zdot - an; %a^n is discrete zdot using bilinear transform
                
                if t > 10000 && drawNR
                    hold on;
                    scatter3(Vrel, z, g1);
                    scatter3(Vrel, z, g2);
                    zlim([-2.5e5, 3.5e5])
                    title("Sample " + t + " Iteration no. " + i);
                    xlabel("$v_{rel}$", 'interpreter', 'latex')
                    ylabel("$z$", 'interpreter', 'latex')
                    set(gca, 'Fontsize', 16)
                    drawnow;
                end
                %% compute derivatives
                
                % Derivative of the steady state function (dz_ss/dv)
                dz_ss = (-2*abs(Vrel) / (strv^2*sig0))*(fs-fc)*espon;
                dz_ssAbs = sign(zss) * dz_ss;
                dalpha_v=0; %d(alpha)/dv
                dalpha_z=0; %d(alpha)/dz
                zss = abs(zss);
                if ((sign(z)==sign(Vrel)) && (abs(z)>z_ba) && (abs(z)<zss))
                    cosarg=cos(sign(z)*arg);
                    dalpha_v=0.5*pi*cosarg*dz_ssAbs*(z_ba-abs(z))/(zss-z_ba)^2;
                    dalpha_z=0.5*pi*cosarg*sign(z)/(zss-z_ba);
                end
                zss = zssPrev;
                
                dzdot_z = -Vrel/zss*(z*dalpha_z + alpha);
                dzdot_v = 1-z * ((alpha+Vrel*dalpha_v)*zss -dz_ss*alpha*Vrel)/zss^2;
                
                %                 dsig1 = 4 * sig0 * (dz_ss * Vrel - zss) / Vrel^2;
                %                 dsig1 = sig1coeff * sig0 * ((dz_ss * alpha * Vrel + zss * (dalpha_v * Vrel - alpha)) / Vrel^2);
                % derivatives of the functions
                %                 dg1v = 2/k + 2 * s0 + (dsig1 * zdot + sig1 * dzdot_v) / (rho * A * h) + sig2 / (rho * A * h);
                dg1v = 2/k + 2 * s0 + sig1 / (rho * A * h) + sig2 / (rho * A * h);
                dg1z = sig0 / (rho * A * h) + sig1 / (rho * A * h) * dzdot_z;
                dg2v = dzdot_v;
                dg2z = dzdot_z - 2/k;
                
                % create Jacobian matrix
                Jac = [dg1v, dg1z; ...
                    dg2v, dg2z];
                determ = dg1v * dg2z - dg1z * dg2v;
                
                % perform vector NR
                prevSolut = [Vrel; z];
                solut = [Vrel; z] - Jac \ [g1; g2];
                VrelCheck = Vrel - 1/determ * (dg2z * g1 - dg1z * g2);
                zCheck = z - 1/determ * (-dg2v * g1 + dg1v * g2);
                VrelPrevIt = Vrel;
                zPrevIt = z;
                Vrel = solut(1);
                z = solut(2);
                
                eps = norm(solut - prevSolut);
                %                 epsTest = 1/2 * ((Vrel - prevSolut(1)) / prevSolut(1) + (z - prevSolut(2)) / prevSolut(2));
                if saveAlphaItFlag
                    alphaIt(i+1) = alpha;
                    vrelIt(i+1) = Vrel;
                    zIt(i+1) = z;
                    zDotIt(i+1) = zdot;
                    epsSave(i+1) = eps;
                end
                %
                %                 scatter(Vrel, z)
                %                 xlim([-2, 2])
                %                 ylim([-1e-3, 1e-3])
                %                 drawnow
                i = i + 1;
                
            end
            %             vRel = vRelTemp;
            %             z = zTemp;
            
            espon=exp(-(Vrel/strv).^2);         %exponential function
            zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
            if Vrel==0
                zss=fs/sig0;
            end
            zssSave(t) = zss;
            
            zssPrev = zss;
            zss = abs(zss);
            alpha = 0;
            
            % elasto-plastic function \alpha (v,z)
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
            
            an = 2/k * (z - zPrev) - anPrev;
            
            %                 zdotPrevIt = zdot;
            
            zdot = Vrel * (1 - alpha * z / zss);
            iSave(t) = i;
            %             sig1 = sig1coeff * sig0 * zss * alpha / Vrel;
            if saveAlphaItFlag
                subplot(5, 1, 1);
                plot(alphaIt);
                title("$\alpha$", 'interpreter', 'latex', 'Fontsize', 16)
                subplot(5, 1, 2);
                plot(vrelIt);
                title("$v_{rel}$", 'interpreter', 'latex', 'Fontsize', 16)
                subplot(5, 1, 3);
                plot(zIt);
                title("$z$", 'interpreter', 'latex', 'Fontsize', 16)
                subplot(5, 1, 4);
                plot(zDotIt);
                title("$\dot{z}$", 'interpreter', 'latex', 'Fontsize', 16)
                subplot(5, 1, 5);
                plot(epsSave);
                title("$\epsilon$", 'interpreter', 'latex', 'Fontsize', 16)
                drawnow;
            end
            
            %             if sig1 * Vrel * alpha / (4 * zss) > sig0
            %                     disp("Not stable")
            %             end
            %             sig1Save(t) = sig1;
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
    uNext = B * u + C * uPrev - excitation / (rho * A / k^2 + s0 / k);
    if bowModel ~= "cos"
        vrelDiff(t) = Vrel - ((uNext(bP) - uPrev(bP)) / (2*k) - Vb);
        bristleEnergy(t) = sig0 / 2 * z^2;
    end
%     Vrel = ((uNext(bP) - uPrev(bP)) / (2*k) - Vb);
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
%     rOCtotBowEnergy(t) = -F * (uNext(bP) - uPrev(bP)) / (2*k);
%     rOCenergy(t) = rOCkinEnergy(t) - rOCpotEnergy(t) - rOCdamp0StringEnergy(t) - rOCdamp1StringEnergy(t) - rOCtotBowEnergy(t);
    if bowModel ~= "cos"
        rOCbowEnergyInput(t) = F * Vb;
        rOCbristleEnergy(t) = sig0 * z * zdot;
        rOCdispEnergy(t) = sig1 * (zdot + 1/2 * Vrel * alpha * z / zss)^2 + z^2 * Vrel * alpha / zss * (sig0 - sig1 * Vrel * alpha /(4*zss)) + sig2 * Vrel^2;
    end
    if bowModel == "elastoPlastic"
        rOCenergy(t) = rOCkinEnergy(t) - rOCpotEnergy(t) - rOCdamp0StringEnergy(t) - rOCdamp1StringEnergy(t) + rOCbowEnergyInput(t) + rOCbristleEnergy(t) + rOCdispEnergy(t) + rOCbowEnergy(t);
    else
        rOCenergy(t) = rOCkinEnergy(t) - rOCpotEnergy(t) - rOCdamp0StringEnergy(t) - rOCdamp1StringEnergy(t);
    end
%         rOCenergy(t) = rOCkinEnergy(t) - rOCpotEnergy(t) - rOCdamp0StringEnergy(t) - rOCdamp1StringEnergy(t) - rOCbowStringEnergy(t);
    
    %% Drawing functions
    if drawString == true && mod(t,drawspeed) == 0 && t >= drawStart
        if drawEnergy == false
            clf
            subplot(4,1,1)
            plot(uNext);
            hold on;
            %         scatter(repmat(floor(bp*N), 20, 1), [-1e-5:1e-6:1e-5-1e-6]+mod(bowVertPos/20,1e-6),'.');
            text(floor(bp*N) - 2, 0, "$V_B =$ " + num2str(Vb, 2), 'interpreter', 'latex', 'horizontalAlignment', 'right');
            title("String displacement at sample " + num2str(t), 'interpreter', 'latex','Fontsize', 16)
            subplot(4,1,2);
            plot(iSave(1:t))
            %             if t > drawStart + 100*drawspeed + 1
            %                 plot(out3(t-100*drawspeed:t-1));
            %             end
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
                %                 alphaIdx = find(round(zVec*1e7) == floor(1e7*z));
                %                 if length(alphaIdx) == 1
                scatter(z, alpha);
                %                     text(floor(1e5*z)*1e-5+3e-5, alphaPlot(alphaIdx) - 0.05, '$z$', 'interpreter', 'latex', 'horizontalAlignment', 'center');
                %                 end
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
                %                 figure('Renderer', 'painters', 'Position', [100 100 800 350])
                subplot(1,1,1)
                numDots = 500;
                mat = [[1:-1/numDots:1/numDots]', [1:-1/numDots:1/numDots]', [1:-1/numDots:1/numDots]'];
                cla;
                if t > numDots
                    plot(vRelSave(t-numDots+1:t), Fsave(t-numDots+1:t), 'k', 'Linewidth', 2);
                    hold on;
                    scatter(vRelSave(t-numDots+1:t), Fsave(t-numDots+1:t), 40, mat);
                end
                ylabel("$f(v,z)$", 'interpreter', 'latex')
                xlabel("$v$", 'interpreter', 'latex')
                title('Force function $f(v,z)$ against relative velocity', 'interpreter', 'latex', 'Fontsize', 16)
                set(gca, 'Fontsize', 20)
                set(gca, 'Position', [0.11 0.15 0.86 0.77])
%                 xlim([-0.42 0.05])
%                 ylim([-2 -1.65])
                grid on;
            end
            drawnow;
            
        elseif drawEnergy == true && t > drawStart - 1000 && mod(t, drawspeed) == 0
            subplot(4,1,1);
            plot(u);
            title("String (state)")
            subplot(4,1,2);
            %             plot(totEnergy(10:t) / totEnergy(10) - 1);
            %             title("Total energy")
            cla;
            if bowModel == "cos"
                plot(totEnergy(10:t)/totEnergy(10) - 1)
            else
                plot(rOCkinEnergy(10:t) - rOCpotEnergy(10:t) - rOCdamp0StringEnergy(10:t) - rOCdamp1StringEnergy(10:t));
                hold on;
                plot(rOCdispEnergy(10:t));
                plot(rOCbowEnergyInput(10:t));
                plot(rOCdispEnergy(10:t));
                plot(rOCbristleEnergy(10:t));
            end
            subplot(4,1,3);
            cla
            plot(rOCenergy(10:t));
            title("Rate of change of energy + bowing energy (should be 0)");
            if bowModel ~= "cos"
                subplot(4,1,4)
                plot(-vrelDiff(10:t));
            end
            drawnow;
        end
    end
    
    % Save output
    out1(t) = uNext(pickup1);
    out2(t) = uNext(pickup2);
    out3(t) = uNext(bP);
    %     if t == 512
    %
    %     end
    % update state vectors
    uPrev = u;
    u = uNext;
end
% figure('Renderer', 'painters', 'Position', [100 100 800 350])
zoom = false;
if zoom
    startPos = Fs / 4;
    endPos = Fs / 4 + Fs / 100;
else
    startPos = 1;
    endPos = lengthSound;
end
subplot(2,1,1)
if zoom
    plot(startPos:endPos, out1(startPos:endPos), 'k', 'Linewidth', 2);
else
    plot(startPos:endPos, out1(startPos:endPos), 'k');
end
title("Output position = 2/3 N")
set(gca, 'Fontsize', 16, 'Position', [0.1 0.59 0.88 0.35])
xlim([startPos, endPos]);
xticks([]);
% yticks([]);
subplot(2,1,2)
if zoom
    plot((startPos:endPos) / Fs, out3(startPos:endPos), 'k', 'Linewidth', 2);
else
    plot((startPos:endPos) / Fs, out3(startPos:endPos), 'k');
end
title("Output position = bowing position (1/4 N)")
set(gca, 'Fontsize', 16, 'Position', [0.1 0.16 0.88 0.35])
xlim([startPos/Fs, endPos/Fs]);
xlabel("Time (s)")

% figure;
% start = 1; endpoint = lengthSound;
% subplot(3,2,1); plot(start:endpoint, out1(start:endpoint))
% title("Output", 'interpreter', 'latex')
% subplot(3,2,2); plot(start:endpoint, zssSave(start:endpoint))
% title("$z_{ss}$", 'interpreter', 'latex')
% subplot(3,2,3); plot(start:endpoint, alphaSave(start:endpoint))
% title("$\alpha(v^n, z^n)$", 'interpreter', 'latex')
% subplot(3,2,4); plot(start:endpoint, zSave(start:endpoint))
% title("$z$", 'interpreter', 'latex')
% subplot(3,2,5); plot(start:endpoint, vRelSave(start:endpoint))
% title("$v_{rel}$", 'interpreter', 'latex')
% subplot(3,2,6);plot(start:endpoint, iSave(start:endpoint))
% title("Number of iterations", 'interpreter', 'latex')
