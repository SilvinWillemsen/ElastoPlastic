%{
%%%%%%%%%%%%%%%%%%%%%%
CONCLUSION
With negative bowing velocity there is no difference between the original
model and when adding absolute values for z_ss
%%%%%%%%%%%%%%%%%%%%%%
%}

clear all;
close all;
clc;

%% Drawing functions
drawString = true;
drawStart = 0;
drawspeed = 10000;
drawEnergy = false;
hysteresis = true;

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
u1Prev = zeros(N,1);
u1 = zeros(N,1);
u1Next = zeros(N,1);

u2Prev = zeros(N,1);
u2 = zeros(N,1);
u2Next = zeros(N,1);

%% Raised cosine
width = 6;
loc = 0.25;
startIdx = floor(floor(loc * N) - width / 2);
endIdx = floor(floor(loc * N) + width / 2);
if bowModel == "cos"
    u1(startIdx : endIdx) = (1 - cos(2 * pi * [0:width] / width)) / 2;
end
u1Prev = u1;

lengthSound = Fs; % Set the length of the output
out1 = zeros(lengthSound, 1);

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

pickup = floor(N/3);    % Pickup position

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
    sig3 = 0.0;
    
    w = (2 * rand(lengthSound,1) - 1);
    
    sigma=[sig0,sig1,sig2];
    
    z1_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!)
    z2_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!)
    
    Vb1Init = -0.1;
    z1Prev = 0;
    z1 = 0;
    
    Vb2Init = -0.1;
    z2Prev = 0;
    z2 = 0;
    tol = 1e-7;
    
    z1dot = Vb1Init;
    z1DotPrev = Vb1Init;
    a1n = z1dot;
    a1nPrev = z1DotPrev;
    
    V1rel = -Vb1Init;
    
    z2dot = Vb2Init;
    z2DotPrev = Vb2Init;
    a2n = z2dot;
    a2nPrev = z2DotPrev;
    
    V2rel = -Vb2Init;
    
    zVec = -z1_ba*5:1e-7:z1_ba*5;
    vRelVec = -0.5:0.001:0.5;
    for ii = 1 : length(vRelVec)
        espon=exp(-(vRelVec(ii)/strv).^2);%exponential function
        zssVec(ii) =sign(vRelVec(ii)).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
        if vRelVec(ii)==0
            zssVec(ii)=fs/sig0;
        end
    end
    z1Save = zeros(lengthSound, 1);
    vRelSave = zeros(lengthSound, 1);
    z1DotSave = zeros(lengthSound, 1);
    F1save = zeros(lengthSound, 1);
    iSave = zeros(lengthSound, 1);
    
else
    a = 100;                % free parameter
    Vb1 = 0.1;               % Bowing speed
    Fb = 1;                % Bowing force / total mass of string
    BM = sqrt(2*a)*exp(1/2);
    tol = 1e-4;
    
    qSave = zeros (lengthSound, 1);
    qPrev = 0;
end

excitation1 = 0;

bowVertPosPrev = 0; % for drawing

scalar = 0;


% w = zeros(lengthSound,1);
% wNoise = awgn(w,10);
% plot(wNoise);
zssVecPlotVal = 1;
% figure('Renderer', 'painters', 'Position', [100 100 900 500])
eVec = 2:N-1;
vec = 3:N-2;
testFlag = false;
for t = 1 : lengthSound
    %     Fn =  t / lengthSound * FnInit
    if bowModel == "elastoPlastic"
        fc=mud*Fn;             % coulomb force
        fs=mus*Fn;             % stiction force
        z1_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!)
    end
%     if t < ramp
%         Vb = VbInit * t / ramp; %* cos(10 * pi * t / Fs);
%     else
        Vb1 = Vb1Init;% * cos(10 * pi * t / Fs);
        Vb2 = Vb2Init;% * cos(10 * pi * t / Fs);
%     end
    % for drawing bow
    %     if drawString == true && bowModel == "elastoPlastic"
    %         bowVertPos = bowVertPosPrev + k/2 * (Vb + VbPrev);
    %         VbPrev = Vb;
    %         bowVertPosPrev = bowVertPos;
    %     end
    if bowModel == "simple"
        b1 = 2/k * Vb1 + 2 * s0 * Vb1 + I' * bB * u1 + I' * bC * u1Prev;
        eps = 1;
        i = 0;
        while eps>tol
            q=qPrev-(1/h*Fb*BM*qPrev*exp(-a*qPrev^2)+2*qPrev/k+2*s0*qPrev+b1)/...
                (1/h * Fb*BM*(1-2*a*qPrev^2)*exp(-a*qPrev^2)+2/k+2*s0);
            eps = abs(q-qPrev);
            qPrev = q;
            i = i + 1;
            if i > 10000
                disp('Nope');
            end
        end
        qSave(t) = q;
        F1save(t) = Fb * sqrt(2 * a) * q * exp(-a*q^2 + 1/2);
        excitation1 = J * F1save(t);
        vRelSave(t) = q;
    elseif bowModel == "elastoPlastic"
        
        % calculate pre-computable part of the FDS
        b1 = (2/k * Vb1 + 2 * s0 * Vb1 + I' * bB * u1 + I' * bC * u1Prev);
        b2 = (2/k * Vb2 + 2 * s0 * Vb2 + I' * bB * u2 + I' * bC * u2Prev);
        if fc>0
            eps = 1;
            i = 0;
            while eps>tol && i < 50
                %% zss and alpha for string 1
                espon=exp(-(V1rel/strv).^2);         %exponential function
                zss1=sign(V1rel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
                if V1rel==0
                    zss1=fs/sig0;
                end
                zss1Save(t) = zss1;
                % elasto-plastic function \alpha (v,z)
                zss1Prev = zss1;
                zss1 = abs(zss1);
                alpha1 = 0;

                if (sign(z1)==sign(V1rel))
                    if ((abs(z1)>z1_ba) && (abs(z1)<zss1))
                        arg1=pi*(z1-sign(z1)*0.5*(zss1+z1_ba))/(zss1-z1_ba);
                        sinarg1=sin(sign(z1)*arg1);
                        alpha1=0.5*(1+sinarg1);
                    elseif (abs(z1)>=zss1)
                        alpha1=1;
                    end
                end
                alpha1Save(t) = alpha1;
                zss1 = zss1Prev;
                z1dotPrevIt = z1dot;
                z1dot = V1rel * (1 - alpha1 * z1 / zss1);
                
                a1n = 2/k * (z1 - z1Prev) - a1nPrev;
                
                %% zss and alpha for string 2
                espon=exp(-(V2rel/strv).^2);         %exponential function
                zss2=sign(V2rel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
                if V2rel==0
                    zss2=fs/sig0;
                end
                zss2Save(t) = zss2;
                
                % elasto-plastic function \alpha (v,z)
                alpha2 = 0;
                
                if (sign(z2)==sign(V2rel))
                    if ((abs(z2)>z2_ba) && (abs(z2)<zss2))
                        arg2=pi*(z2-0.5*(zss2+z2_ba))/(zss2-z2_ba);
                        sinarg2=sin(arg2);
                        alpha2=0.5*(1+sinarg2);
                    elseif (abs(z2)>=zss2)
                        alpha2=1;
                    end
                end
                alpha2Save(t) = alpha2;
                z2dotPrevIt = z2dot;
                z2dot = V2rel * (1 - alpha2 * z2 / zss2);
                
                a2n = 2/k * (z2 - z2Prev) - a2nPrev;
                
                
                %% COMPARE
%                 if t>291
%                     zss1 + zss2
% %                     alpha1-alpha2
% %                     z1+z2
% %                     V1rel + V2rel
% %                     testFlag = true;
%                 end
                
                %% compute derivatives for string 1
                
                % Derivative of the steady state function (dz_ss/dv)
                d1z_ss = (-2*V1rel*sign(V1rel) / (strv^2*sig0))*(fs-fc)*espon;
                d1z_ssAbs = sign(zss1) * d1z_ss;
                d1alpha_v=0; %d(alpha)/dv
                d1alpha_z=0; %d(alpha)/dz
                zss1 = abs(zss1);
                if ((sign(z1)==sign(V1rel)) && (abs(z1)>z1_ba) && (abs(z1)<zss1))
                    cosarg1=cos(sign(z1)*arg1);
                    d1alpha_v=0.5*pi*cosarg1*sign(zss1)*d1z_ssAbs*(z1_ba-sign(z1)*z1)/(zss1-z1_ba)^2;
                    d1alpha_z=0.5*pi*cosarg1*sign(z1)/(zss1-z1_ba);
                end
                zss1 = zss1Prev;
                
                %                 zdot = Vrel * (1-alpha * z / zss); %non-linear function estimate
                d1fnl_z = -V1rel/zss1*(z1*d1alpha_z + alpha1);
                d1fnl_v = 1-z1 * ((alpha1+V1rel*d1alpha_v)*zss1-d1z_ss*alpha1*V1rel)/zss1^2;
                
                % functions to perform newton raphson on
                g11 = (2/k + 2 * s0 + sig2 / (rho * A * h)) * V1rel + (sig0 * z1 + sig1 * a1n ...
                    + sig3 * w(t)) / (rho * A * h) + b1;
                g12 = z1dot - a1n; %a^n (discrete zdot, bilinear transform)
                
                % derivatives of the functions
                d1g1v = 2/k + 2 * s0 + sig2 / (rho * A * h);
                d1g1z = sig0 / (rho * A * h) + sig1 / (rho * A * h) * 2/k;
                d1g2v = d1fnl_v;
                d1g2z = d1fnl_z - 2/k;
                
                % create Jacobian matrix
                Jac1 = [d1g1v, d1g1z; d1g2v, d1g2z];
                
                %% compute derivatives for string 2
                
                % Derivative of the steady state function (dz_ss/dv)
                d2z_ss = (-2*V2rel*sign(V2rel) / (strv^2*sig0))*(fs-fc)*espon;
                d2alpha_v=0; %d(alpha)/dv
                d2alpha_z=0; %d(alpha)/dz
                if ((sign(z2)==sign(V2rel)) && (abs(z2)>z2_ba) && (abs(z2)<zss2))
                    cosarg2=cos(arg2);
                    d2alpha_v=0.5*pi*cosarg2*d2z_ss*(z2_ba-z2)/(zss2-z2_ba)^2;
                    d2alpha_z=0.5*pi*cosarg2/(zss2-z2_ba);
                end
                
                %                 zdot = Vrel * (1-alpha * z / zss); %non-linear function estimate
                d2fnl_z = -V2rel/zss2*(z2*d2alpha_z + alpha2);
                d2fnl_v = 1-z2 * ((alpha2+V2rel*d2alpha_v)*zss2-d2z_ss*alpha2*V2rel)/zss2^2;
                
                % functions to perform newton raphson on
                g21 = (2/k + 2 * s0 + sig2 / (rho * A * h)) * V2rel + (sig0 * z2 + sig1 * a2n ...
                    + sig3 * w(t)) / (rho * A * h) + b2;
                g22 = z2dot - a2n; %a^n (discrete zdot, bilinear transform)
                
                % derivatives of the functions
                d2g1v = 2/k + 2 * s0 + sig2 / (rho * A * h);
                d2g1z = sig0 / (rho * A * h) + sig1 / (rho * A * h) * 2/k;
                d2g2v = d2fnl_v;
                d2g2z = d2fnl_z - 2/k;
                
                % create Jacobian matrix
                Jac2 = [d2g1v, d2g1z; d2g2v, d2g2z];
                
                
                %% perform vector NR string 1
                solut1 = [V1rel; z1] - Jac1 \ [g11; g12];
                V1relSave(t) = V1rel;
                V1rel = solut1(1);
                z1 = solut1(2);

                
                %% perform vector NR string 2
                solut2 = [V2rel; z2] - Jac2 \ [g21; g22];
                V2relSave(t) = V2rel;
                V2rel = solut2(1);
                z2 = solut2(2);

                test = V1rel * (1 - alpha1 * z1 / zss1);
                eps = abs(test-z1dot);
                
                i = i + 1;
                
            end
            iSave(t) = i;
            F1 = (sig0 * z1 + sig1 * a1n + sig2 * V1rel + sig3 * w(t));
            F1save(t) = F1;
            F2 = (sig0 * z2 + sig1 * a2n + sig2 * V2rel + sig3 * w(t));
            F2save(t) = F2;
        else
            z1dot=0; err=0; count=0;
            z1=0; f_fr=0; v=vs;
        end
        excitation1 = J*F1;
        VRel1SaveNow(t) = V1rel;
        z1Save(t) = z1;
        z1DotSave(t) = z1dot;
%                         plot(VrelSave(1:t))
%                 drawnow;
        z2DotSave(t) = z2dot;
        excitation2 = J*F2;
        VRel2SaveNow(t) = V2rel;
        z2Save(t) = z2;
        z2DotSave(t) = z2dot;
%         if t == 774
%             load VRelSave.mat
%             plot(VRelSaveNow);
%             hold on;
%             plot(VrelSave);
%         end
        % update z and zdot
        z1Prev = z1;
        z1DotPrev = z1dot;
        a1nPrev = a1n;
        
        z2Prev = z2;
        z2DotPrev = z2dot;
        a2nPrev = a2n;
        
    elseif bowModel == "hyperbolic"
        v = 1/k * (I' * u1 - I' * u1Prev);
        mu = (mud + (mus - mud) * Vb1/2) / (Vb1/2 + v - Vb1);
        excitation1 = J*Fn*mu;
    else
        excitation1 = 0;
    end
    
    %% Update FDS
    u1Next = (B * u1 + C * u1Prev) - excitation1 / (rho * A / k^2 + s0 / k);
    u2Next = (B * u2 + C * u2Prev) - excitation2 / (rho * A / k^2 + s0 / k);

    %% Calculate energy of the string
    kinEnergy(t) = rho * A / 2 * h * sum((1/k * (u1 - u1Prev)).^2);
    potEnergy(t) = T / 2 * 1/h * sum((u1(3:N) - u1(2:N-1)) .* (u1Prev(3:N) - u1Prev(2:N-1)))...
        + E * Iner / 2 * 1/h^3 * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1)) ...
        .* (u1Prev(eVec+1) - 2 * u1Prev(eVec) + u1Prev(eVec-1)));
    totEnergy(t) = kinEnergy(t) + potEnergy(t);
    
    
    %% Calculate Rate of Changes of Energy
    rOCkinEnergy(t) = h * rho * A / (2 * k^3) * sum((u1Next - 2 * u1 + u1Prev) .* (u1Next - u1Prev));
    rOCpotEnergy(t) = h * T / (2*k*h^2) * sum((u1(vec+1) - 2 * u1(vec) + u1(vec-1)).* (u1Next(vec) - u1Prev(vec))) ...
        - h * E * Iner / (2 * k * h^4) * sum((u1(vec+2) - 4 * u1(vec+1) + 6 * u1(vec) - 4 * u1(vec-1) + u1(vec-2)) .* (u1Next(vec) - u1Prev(vec)));%...
    rOCdamp0StringEnergy(t) = -2 * s0 * h / (4 * k^2) * sum((u1Next - u1Prev).*(u1Next - u1Prev));
    rOCdamp1StringEnergy(t) = 2 * h * s1 / (2 * k^2 * h^2) * sum((u1(eVec+1) - 2 * u1(eVec) + u1(eVec-1) - u1Prev(eVec+1) + 2 * u1Prev(eVec) - u1Prev(eVec-1)) .* (u1Next(eVec) - u1Prev(eVec)));
    rOCbowStringEnergy(t) = -h * sum(I .* excitation1) * (u1Next(bP) - u1Prev(bP)) / (2 * k);
    rOCenergy(t) = rOCkinEnergy(t) - rOCpotEnergy(t) - rOCdamp0StringEnergy(t) - rOCdamp1StringEnergy(t) - rOCbowStringEnergy(t);
    
    %% Drawing functions
    if drawString == true && mod(t,drawspeed) == 0 && t > drawStart
        if drawEnergy == false 
            clf
            title("Sample " + t)
            subplot(4,1,1)
            plot(u1Next);
            hold on;
            plot(u2Next);
            %         scatter(repmat(floor(bp*N), 20, 1), [-1e-5:1e-6:1e-5-1e-6]+mod(bowVertPos/20,1e-6),'.');
            text(floor(bp*N) - 2, 0, "$V_B =$ " + num2str(Vb1, 2), 'interpreter', 'latex', 'horizontalAlignment', 'right');
            title("String displacement at sample " + num2str(t), 'interpreter', 'latex','Fontsize', 16)
            if hysteresis == false
                %% Steady State Curve
                subplot(2,4,5)
                plot(vRelVec, zssVec)
                hold on;
                zssVecVal = find (round(vRelVec * 1e3) == round(V1rel * 1e3));
                if length(zssVecVal) == 1
                    zssVecPlotVal = zssVecVal;
                end
                scatter(round(V1rel * 1e3)*1e-3, zssVec(zssVecPlotVal));
                text(round(V1rel * 1e3)*1e-3, zssVec(zssVecPlotVal) - 1.5e-4, "$v =$ " + num2str(V1rel, 2), 'interpreter', 'latex', 'horizontalAlignment', 'center');
                xlim([-0.5 0.5])
                %         title("$z$", 'interpreter', 'latex')
                xlabel('$v$','interpreter', 'latex')
                ylabel("$z_{ss}(v)$", 'interpreter', 'latex')
                title('Steady-state curve $z_{ss}(v)$', 'interpreter', 'latex', 'Fontsize', 16)
                
                %% Adhesion Map
                subplot(2,4,6)
                alphaPlot = zeros(length(zVec),1);
                zssPlot = abs(zss1);
                for ii = 1:length(zVec)
                    if ((abs(zVec(ii))>z1_ba) && (abs(zVec(ii))<zssPlot))
                        arg1=pi*(zVec(ii)-sign(zVec(ii))*0.5*(zssPlot+z1_ba))/(zssPlot-z1_ba);
                        sinarg1=sin(sign(zVec(ii))*arg1);
                        alphaPlot(ii)=0.5*(1+sinarg1);
                    elseif (abs(zVec(ii))>=zssPlot)
                        alphaPlot(ii)=1;
                    end
                end
                plot(zVec, alphaPlot);
                text(zss1 + 5e-5, 0.5, '$z_{ss}(v)$', 'interpreter', 'latex', 'horizontalAlignment', 'center')
                
                hold on;
                plot([zss1 zss1], [min(ylim) max(ylim)], '--');
                alphaIdx = find(round(zVec*1e7) == floor(1e7*z1));
                if length(alphaIdx) == 1
                    scatter(floor(1e5*z1)*1e-5, alphaPlot(alphaIdx));
                    text(floor(1e5*z1)*1e-5+3e-5, alphaPlot(alphaIdx) - 0.05, '$z$', 'interpreter', 'latex', 'horizontalAlignment', 'center');
                end
                %         title('$\alpha$','interpreter', 'latex', 'Fontsize', 18)
                xlabel('$z$','interpreter', 'latex')
                ylabel('$\alpha(v,z)$','interpreter', 'latex')
                title('Adhesion map $\alpha(v,z)$', 'interpreter', 'latex', 'Fontsize', 16)
                
                %% Velocity of Mean Bristle Displacement
                subplot(2,4,7)
                plot(z1DotSave(1:t))
                xlabel('$n$ (samples)','interpreter', 'latex')
                ylabel("$\dot z$", 'interpreter', 'latex')
                title('Velocity of the mean bristle displacement', 'interpreter', 'latex', 'Fontsize', 16)
                
                %% Mean Bristle Displacement
                subplot(2,4,8)
                plot(z1Save(1:t))
                title("$z$", 'interpreter', 'latex')
                xlabel('$n$ (samples)','interpreter', 'latex')
                ylabel("$z$", 'interpreter', 'latex')
                title('Mean bristle displacement $z$', 'interpreter', 'latex', 'Fontsize', 16)
            else 
                subplot(4,1,2)
                plot(u1Next+u2Next);
                subplot(4,1,3)
%                 if t > drawspeed * 3
%                     cla
                    plot(alpha1Save(1:t))
                    hold on;
                    plot(alpha2Save(1:t))
%                 end
                subplot(4,1,4)
                plot(alpha1Save(1:t) - alpha2Save(1:t))
%                 numDots = 500;
%                 mat = [[1:-1/numDots:1/numDots]', [1:-1/numDots:1/numDots]', [1:-1/numDots:1/numDots]'];
%                 %             vec = t-numDots + 1 : t;
%                 cla;
%                 if t > numDots
%                                    plot(vRelSave(t-numDots+1:t), F1save(t-numDots+1:t), 'k', 'Linewidth', 2);
%                                hold on;
%                     scatter(vRelSave(t-numDots+1:t), F1save(t-numDots+1:t), 40, mat);
%                 end
%                 ylabel("$f(v,z)$", 'interpreter', 'latex')
%                 xlabel("$v$", 'interpreter', 'latex')
%                 title('Force function $f(v,z)$ against relative velocity', 'interpreter', 'latex', 'Fontsize', 16)
%                 %             xlim([-0.25 0.05])
%                 %             ylim([-0.2 1])
%                 set(gca, 'Fontsize', 20)
%                 %             set(gca, 'Position', [0.08 0.11 0.89 0.82])
%                 grid on;
            end
            drawnow;
            
        elseif drawEnergy == true && t > drawStart - 1000 && mod(t, drawspeed) == 0
            subplot(3,1,1);
            plot(u1);
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
    out1(t) = u1Next(pickup);
    out2(t) = u2Next(pickup);
    %     out2(t) = u(pickup) + 1 / 2 * (uNext(pickup) - uPrev(pickup));
    
    % update state vectors
    u1Prev = u1;
    u1 = u1Next;
    u2Prev = u2;
    u2 = u2Next;
end

plot(out1);
hold on;
plot(out2);
figure;
plot(out1 + out2')