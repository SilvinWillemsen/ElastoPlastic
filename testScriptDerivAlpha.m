%{
Initialise the variables using the StiffStringElastoPlastic.m file and
stop right before the time-loop.
%}


arg = 0;
Vrel = [ones(500,1) * -0.1; ones(501, 1)*0.1];
zVec = -z_ba*5:z_ba/100:z_ba*5;

espon=exp(-(Vrel/strv).^2);         %exponential function
zss=sign(Vrel).*(fc +(fs-fc)*espon)/sig0;   %steady state curve: z_ss(v)
if Vrel==0
    zss=fs/sig0;
end


alpha = zeros(length(zVec), 1);
dalpha_v=zeros(length(zVec), 1); %d(alpha)/dv
dalpha_z=zeros(length(zVec), 1); %d(alpha)/dz

zssPrev = zss;
zss = abs(zss);
for j = 1:length(zVec)
    if (sign(zVec(j))==sign(Vrel(j)))
        if ((abs(zVec(j))>z_ba) && (abs(zVec(j))<zss(j)))
            arg=pi*(zVec(j)-sign(zVec(j))*0.5*(zss(j)+z_ba))/(zss(j)-z_ba);
            sinarg=sin(sign(zVec(j))*arg);
            in = 0.5*(1+sinarg);
            alpha(j)=in;
        elseif (abs(zVec(j))>=zss(j))
            alpha(j)=1;
        end
    end
    
    dz_ss = (-2*abs(Vrel(j)) / (strv^2*sig0))*(fs-fc)*espon;
    dz_ssAbs = sign(zss(j)) * dz_ss(j);
    zss(j) = abs(zss(j));
    if ((sign(zVec(j))==sign(Vrel(j))) && (abs(zVec(j))>z_ba) && (abs(zVec(j))<zss(j)))
        cosarg=cos(sign(zVec(j))*arg);
        dalpha_v(j)=0.5*pi*cosarg*dz_ssAbs*(z_ba-abs(zVec(j)))/(zss(j)-z_ba)^2;
        dalpha_z(j)=0.5*pi*cosarg*sign(zVec(j))/(zss(j)-z_ba);
    end
end

subplot(3, 1, 1);
plot(alpha)
title("\alpha")
subplot(3, 1, 2);
plot(dalpha_v)
title("$\frac{\partial\alpha}{\partial v}$", 'interpreter', 'latex');
subplot(3, 1, 3);
plot(dalpha_z)
title("$\frac{\partial\alpha}{\partial z}$", 'interpreter', 'latex');