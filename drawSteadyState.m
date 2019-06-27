clear all;
close all;
%%%% the Contact Force (be with you) %%%%%%%%%
mus=0.4;               % static friction coeff
mud=0.2;              % dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
strv=0.1;             % "stribeck" velocity

plotzss = true;

Fn = 5;

fc=mud*Fn;             % coulomb force
fs=mus*Fn;             % stiction force
ssparam=[fc,fs,strv];            

sig0=1e4;           % bristle stiffness
sig1=.1*sqrt(sig0);   % bristle damping
sig2=0.4;            % viscous friction term 
sigma=[sig0,sig1,sig2];

z_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!) 


if plotzss
    Vrel = -0.5:0.001:0.5;
else
    Vrel = [0.2 * ones(3001,1)];
end
espon=exp(-(Vrel/strv).^2);         %exponential function
zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
if Vrel==0
  zss=fs/sig0;
end

if plotzss
    figure('Renderer', 'painters', 'Position', [100 100 600 250])
    set(gca, 'Position', [0.08 0.18 0.91 0.73])
    plot(Vrel, zss, 'k', 'LineWidth', 2)
    title("Steady-state function")
    xlabel('$v$', 'interpreter', 'latex')
    ylabel('$z_{ss}(v)$', 'interpreter', 'latex')
    ylim([-2.5e-4 2.5e-4])
    grid on;
    set(gca, 'Fontsize', 16, 'Linewidth', 2);
    return;
end

z = -0.00015:0.0000001:0.00015;
alphaVar=zeros(3001,1);
zssPrev = zss;
zss = abs(zss);
for i = 1:3001
%     if (sign(z(i))==sign(Vrel(i)))
        if ((abs(z(i))>z_ba) && (abs(z(i))<zss(i)))
            arg=pi*(z(i) - 0.5*(zss(i)+z_ba))/(zss(i)-z_ba);
            sinarg(i)=sin(arg);
            alphaVar(i)=0.5*(1+sinarg(i));
        elseif (abs(z(i))>zss(i))
            alphaVar(i)=1;
        end
%     end
end
zss = zssPrev;

% numSamps = length(z);
% alphaVar=zeros(numSamps,1);
% for i = 1:numSamps
%     if ((abs(z(i))>z_ba) && (abs(z(i))<zss))
%         arg=pi*(z(i)-sign(z(i)) * 0.5*(zss+z_ba))/(zss-z_ba);
%         sinarg=sin(sign(z(i)) * arg);
%         alphaVar(i)=0.5*(1+sinarg);
%     elseif (abs(z(i))>zss)
%         alphaVar(i)=1;
%     end
% end

figure('Renderer', 'painters', 'Position', [100 100 600 420])
hold on;
ylim([-0.07 1.07])
xlim([-1.5e-4 1.5e-4])
patch([-z_ba z_ba z_ba -z_ba], [max(ylim) max(ylim) min(ylim) min(ylim)], [1 1 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch([z_ba zss(1) zss(1) z_ba], [max(ylim) max(ylim) min(ylim) min(ylim)], [1 0.65 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch([zss(1) zss(1)*1.5 zss(1)*1.5 zss(1) ], [max(ylim) max(ylim) min(ylim) min(ylim)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch(-[z_ba zss(1) zss(1) z_ba], [max(ylim) max(ylim) min(ylim) min(ylim)], [1 0.65 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
patch(-[zss(1) zss(1)*1.5 zss(1)*1.5 zss(1) ], [max(ylim) max(ylim) min(ylim) min(ylim)], [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
locations = [-zss(1); -z_ba; 0; z_ba; zss(1)];
for i = 1:length(locations)
    if i == 3
        plot([locations(i) locations(i)], [min(ylim) max(ylim)], 'Linewidth', 2, 'color', [0.6 0.6 0.6], 'Linestyle', '--');
    else
        plot([locations(i) locations(i)], [min(ylim) max(ylim)], 'Linewidth', 2, 'color', [0.6 0.6 0.6]);
    end
end
for i = 0:0.1:1
    plot([min(xlim) max(xlim)], [i i], 'Linewidth', 2, 'color', [0.6 0.6 0.6]);
end
fontsize = 16;
xNames = {'-$z_{ss}(v)$', '-$z_{ba}$', '$0$', '$z_{ba}$', '$z_{ss}(v)$'};
text(locations, repmat(-0.11, 5, 1), xNames, 'horizontalAlignment', 'center', 'Fontsize', 18, 'interpreter', 'latex')
text(0, 0.85, '$|z|\leq z_{ba}$', 'interpreter', 'latex', 'Fontsize', fontsize, 'horizontalAlignment', 'center')
text((z_ba+zss(1)) / 2, 0.85, '$z_{ba}<|z|<|z_{ss}(v)|$', 'interpreter', 'latex', 'Fontsize', fontsize, 'horizontalAlignment', 'center', 'rotation', -90)
text((zss(1) + 1.5e-4) / 2, 0.85, '$|z|\geq|z_{ss}(v)|$', 'interpreter', 'latex', 'Fontsize', fontsize, 'horizontalAlignment', 'center')

set(gca, 'Fontsize', 15, 'TickLabelInterpreter', 'latex')
text(0, -0.15, '$z$','interpreter', 'latex', 'Fontsize', 20, 'horizontalAlignment', 'center')
title('Plot for $\alpha(v,z)$ if sgn($v$) == sgn($z$)','interpreter', 'latex', 'Fontsize', 20)
ylabel('$\alpha(v,z)$','interpreter', 'latex')
plot(z, alphaVar, 'k', 'Linewidth', 2);
xticks([])
%%DRAW VERTICAL GRIDLINES
% scatter([-z_ba; z_ba; -zss(1); zss(1)], [0; 0; 1; 1], 50, 'circle')

set(gca, 'Position', [0.08 0.09 0.91 0.84])

patch([min(xlim) max(xlim) max(xlim) min(xlim)], [max(ylim) max(ylim) min(ylim) min(ylim)], [0 0 0], 'FaceColor', 'none', 'EdgeColor', 'black', 'LineWidth', 2);
% grid on;
set(gca, 'Linewidth', 2, 'GridAlpha', 0.5)
