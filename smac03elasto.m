% bowdsprings: two multi-mass resonators interact through a % non-linear dynamic elasto-plastic friction force (wow :-)% see Dupont et al. "Single state elasto-plastic friction models",% IEEE Trans. Autom. Control, june 2002%% authors: SObbers (alphabetic: F. Avanzini, D. Rocchesso, S. Serafin)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear all;close all;clc;global Fs;%%%%%% compatibility MATLAB/Octave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                 %f_prog=0;             % program flag: 0 MATLAB, 1 Octave%if f_prog==0  myred='''r''';  mygreen='''g''';  myblue='''b''';  %myblack='''k''';  mymagenta='''m''';    mygridon='grid on';  mygridoff='grid off';  myclf='clf';  myreplot='';% else  myred='1';  mygreen='2';  myblue='3';  mymagenta='4';  mygridon='grid ("on")';  mygridoff='grid ("off")';  myclf='clearplot';  myreplot='replot';%end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fs=44100;             % sampling rateh = 2*Fs;             % bilinear transform constantsmplength=h;      % sound length in samples% (for clarity the first object is called "b", bow%   and the second one is called "r", resonator)%%%%%%%%%%%% control parameters %%%%%%%%%%Vin= 0;           % initial bow velocityVb= .5           % steady state bow velocity                    fn=1;             % normal force%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% physical parameters%%%% normalized spectra for various objects %%%%%%%%%%npartials=15;ffbar=1/3.011^2*[3.011^2, 5^2, 7^2, 9^2, 11^2, 13^2, 15^2 17^2, 19^2, 21^2, ...      23^2, 25^2, 27^2, 29^2, 31^2];   % normalized partials for a free-free barcfbar=1/1.194^2*[1.194^2, 2.988, 5^2, 7^2, 9^2, 11^2, 13^2, 15^2 17^2, 19^2, ...      21^2, 23^2, 25^2, 27^2, 29^2];   % normalized partials for a clamped-free barffstring=[1:15];                    % normalized partials for a fixed-fixed string%%%% the Bow %%%%%%%%%%%%%%%%%%%%object_b=ffstring; % type of object (see above)SPRING_b=1;  % no. of partialspitch_b=0; % frequency of the first partial [Hz]qfact_b=1e5;      % quality factor of each 2nd order oscillatormass_b=50e-3;    % mass of each 2nd order oscillator%object_b=ffstring; % type of object (see above)%SPRING_b=3;  % no. of partials%pitch_b=4000; % frequency of the first partial [Hz]%qfact_b=1;      % quality factor of each 2nd order oscillator%mass_b=10e-3;    % mass of each 2nd order oscillatorSPRING_b = min(SPRING_b,npartials);partials_b=pitch_b*[object_b(1:SPRING_b)]';m_b = mass_b*ones(SPRING_b,1);              % oscillator massesq_b = qfact_b*ones(SPRING_b,1);             % quality factorsk_b = m_b.*(2*pi*partials_b).^2;      % oscillator elastic constantsr_b = 2*pi*partials_b.*m_b./q_b;        % oscillator damping coefficientt_eb = 2 * m_b ./ r_b%%%% the Resonator %%%%%%%%%%%%%%%%%%%%object_r=ffstring; % type of object (see above)%SPRING_r=5;     % no. of partials%pitch_r=440;   % frequency of the first partial [Hz]SPRING_r=2;     % no. of partialspitch_r=140;   % frequency of the first partial [Hz]qfact_r=500;      % quality factor of each 2nd order oscillatormass_r=1e-3;     % mass of each 2nd order oscillatorSPRING_r = min(SPRING_r,npartials);partials_r=pitch_r*[object_r(1:SPRING_r)]';m_r = mass_r*ones(SPRING_r,1);              % oscillator massesq_r = qfact_r*ones(SPRING_r,1);             % quality factorsk_r = m_r.*(2*pi*partials_r).^2;      % oscillator elastic constantsr_r = 2*pi*partials_r.*m_r./q_r;        % oscillator damping coefficientst_e = 2 * m_r ./ r_r%%%% the Contact Force (be with you) %%%%%%%%%mus=.4;                %static friction coeffmud=.2;                %dynamic friction coeff (must be < mus!!)strv=1e-1;             % "stribeck" velocityfc=mud*fn;             % coulomb forcefs=mus*fn;             % stiction forcessparam=[fc,fs,strv];sig0=1e4;           % bristle stiffnesssig1=.1*sqrt(sig0);   % bristle dampingsig2=0.4;            % viscous friction term sigma=[sig0,sig1,sig2];z_ba=0.7*fc/sig0;    % break-away displacement (has to be < f_c/sigma_0!!) %fe_b=.5;fe_b=fc +(fs-fc)*exp(-(Vb/strv)^2) +sig2*Vb% this is the external force on the bow such that the bow% steady-state velocity (after an initial transient) is Vbfe_r=0;  %no external forces on the resonator%%%%% initializations %%%%%%%%%%%%%%%%%%%%%%%x_b1=zeros(1,SPRING_b);         % initial bow displacs               dotx_b1=Vin*ones(1,SPRING_b);     % initial bow velocities%x_bs=zeros(1,SPRING_b);         % initial bow displacs (history)%dotx_bs=V*ones(1,SPRING_b);     % initial bow velocities (history)x_r1=zeros(1,SPRING_r);         % initial res. displacsdotx_r1=zeros(1,SPRING_r);      % initial res. velocitiesx_rs=zeros(1,SPRING_r)';         % initial res. displacs (history)dotx_rs=zeros(1,SPRING_r)';      % initial res. velocities (history)%v1=0;               % initial relative velocityz1=0;               % initial mean bristle displacementy1=0;               % initial non-linear functionf_fr1=0;                % initial friction force f_tot_r1=0;         % initial total force on resonatorf_tot_b1=0;         % initial total force on bow%%%% resonator dynamics %%%%% f_tot_r= fe_r +f_fr;      <--- total force acting on each mass of the res.% m \ddotx + r \dotx + k x = f_tot_r;   <--- diff. eq. for each mass of the res.% state variable (matrix) representation% _/  m \dot(dotx) + r dotx + k x = f_tot_r; <--- for each mass of the res.%  \  \dot(x) = dotx;% Replacing \dot with h\frac{1 - z^{-1}}{1 + z^{-1}}:% x(n) = a_r [x(n-1), dotx(n-1), f_tot_r(n-1)]^T + den_r*f_tot_r(n);% dotx(n) = b_r [x(n-1), dotx(n-1), f_tot_r(n-1)]^T + h*den_r*f_tot_r(n);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%den_r = 1./(h^2*m_r + h*r_r + k_r);%printf("den_r (size m_r = %f): h=%f, m_r=%f, r_r=%f, k_r=%f; t_e = 2*m_r/r_r = %f \n", length(m_r), h, m_r(1), r_r(1), k_r(1), 2*m_r(1)/r_r(1));%printf("1/(%f + %f +%f) - - t_e = %f \n", h^2*m_r, h*r_r, k_r, m_r/r_r);a_r=diag(den_r)*[(1./den_r-2*k_r), 2*h*m_r,ones(SPRING_r,1)];b_r=diag(den_r)*h*[-2*k_r,(2*h*m_r-1./(h*den_r)),ones(SPRING_r,1)];%%%% bow dynamics %%%%%%%%% f_tot_b= fe_b -f_fr; <--- total force acting on each mass of the bow% ... same as above (with obvious changes of notation)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%den_b = 1./(h^2*m_b + h*r_b + k_b);%printf("1./(h^2*m_b + h*r_b + k_b) = 1/(%f + %f +%f);  \n", h^2*m_b(1), h*r_b(1), k_b(1));%printf("den_b: h=%f, m_b=%f, r_b=%f, k_b=%f;  \n", h(1), m_b(1), r_b(1), k_b(1));a_b=diag(den_b)*[(1./den_b-2*k_b), 2*h*m_b,ones(SPRING_b,1)];b_b=diag(den_b)*h*[-2*k_b,(2*h*m_b-1./(h*den_b)),ones(SPRING_b,1)];%%%% bristles dynamics %%%%%%%% \dotz = v(1-\alpha(v,z)*z/z_{ss}(v))= y  <--- diff. eq. for mean bristle displ.% Replacing \dot with h\frac{1 - z^{-1}}{1 + z^{-1}}:% z(n) = a_z [z(n-1), y(n-1)]^T + 1/h * y(n)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%a_z=[1,1/h];%%%% nonlinear friction (see nlfriction.m) %%%%%%%%%%%%%den_rden_b%printf("entryK2*1E7 = %f, entryK1*1E7 = %f \n", sum(den_r)*1E7, sum(den_b)*1E7);b=h*sum(den_r) +h*sum(den_b)sig0sig1sig2hK=[-b/(1+sig2*b)*(sig0/h +sig1), 1/h]  % K matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   pause;%%%%%%%%%%%%%%%%%%%% sample-by-sample loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  some signals %%%%%time=(1:smplength)/Fs; xrt=zeros(1,smplength);      % resonator displacement xbt=zeros(1,smplength);      % bow displacementvrt=zeros(1,smplength);      % resonator velocity vbt=zeros(1,smplength);      % bow velocityft=zeros(1,smplength);        % contact forcezt=zeros(1,smplength);      % mean bristle displacementcountt=zeros(1,smplength);      %errt=zeros(1,smplength);      %%=zeros(1,smplength);      %%=zeros(1,smplength);      %%%%%%%%%%% a counter%%%%%%%%%%%fprintf(1,'\n Computation n. %d, %d, %d',ipitch,ik,iq);%pause;percent = 0;m1 = 0;m2 = 0;N100 = round(smplength/100);N10 = round(smplength/10);fprintf(1,'\n percentage computed.... \n  0');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for i = 1:smplength      past_b=[x_b1;dotx_b1;f_tot_b1*ones(1,SPRING_b)];   % vector of past values (bow dynamics)   for ib=1:SPRING_b         x_bs(ib)=a_b(ib,:)*past_b(:,ib);          % computable part      dotx_bs(ib)=b_b(ib,:)*past_b(:,ib);       % of bow dynamics   end      past_r=[x_r1;dotx_r1;f_tot_r1*ones(1,SPRING_r)];  % matrix of past values (res. dynamics)   for ir=1:SPRING_r         x_rs(ir)=a_r(ir,:)*past_r(:,ir);          % computable part of      dotx_rs(ir)=b_r(ir,:)*past_r(:,ir);       % resonator dynamics   end      %%%%% bristles %%%%%%%   past_z=[z1;y1];             % vector of past values (bristle dynamics)   zs=a_z*past_z;              % computable part of                               % bristle dynamics                                  vs= 1/(1 +sig2*b)* ...     %computable part of relative velocity      (sum(dotx_bs) -sum(dotx_rs) -b*sig0*zs +h*sum(den_b)*fe_b -h*sum(den_r)*fe_r);            %%%% nonlinear friction %%%%%%%%%%%%%   if fc>0      [y,err,count] = nlfriction(vs,zs,y1,ssparam,z_ba,sigma,K);      v=vs+K(1)*y;                   % relative velocity      z=zs+K(2)*y;                   % mean bristle displacement      f_fr=sig0*z +sig1*y +sig2*v;   % the friction force!   else                              % negative normal force => no contact      y=0; err=0; count=0;      z=0; f_fr=0; v=vs;   end     f_tot_b=fe_b -f_fr;                   % total force on bow   f_tot_r=fe_r +f_fr;                  % total force on resonator        %%%%% displacements and velocities %%%%%   x_b=x_bs' +den_b'*f_tot_b;   dotx_b=dotx_bs' +h*den_b'*f_tot_b;      x_r = x_rs' + den_r'*f_tot_r;   dotx_r=dotx_rs'+ h*den_r'*f_tot_r;      %%%% state update %%%%%%%%%%%%   y1=y;  f_tot_b1=f_tot_b;  f_tot_r1 = f_tot_r;   x_b1=x_b;  dotx_b1=dotx_b;   x_r1=x_r;  dotx_r1=dotx_r;   z1=z;   %v1=v;         %%%% output signals %%%%%%%%%%   ft(i) = f_fr;   xbt(i) = sum(x_b);    xrt(i) = sum(x_r);   vbt(i) = sum(dotx_b);    vrt(i) = sum(dotx_r);   zt(i)  = z;   countt(i) = count;   errt(i) = err;      %%%%%%%% a counter %%%%%%%%%      m1 = m1+1;      m2 = m2+1;      if (m1>N100 & m2<N10)         m1=0;         fprintf(1,'.');      end      if m2>N10         m1=0;         m2=0;         percent = percent+10;         fprintf(1,'\n %d',percent);      end      %%%%%%%%%%%%%%%%%%%%%%%%%%%%   end%%%%%%%%%%%% write wav file  %%%%%%%%%%%wavname='object.wav';if f_prog==0%  wavwrite(vrt/(1.001*max(abs(vrt))),Fs,wavname);else % ausave(wavname, vrt'/(1.001*max(abs(vrt))), Fs, 'short');end;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%final plotswin = 1:smplength;figure(1);subplot(2,2,1)plot(time(win),xrt(win),'r');hold on;plot(time(win),xbt(win),'g');grid;xlabel('time [s]')ylabel('displacs [m]')axis([0, time(max(size(win))), -0.0015, max(xbt)+0.0015]);hold off;myreplot;subplot(2,2,2)plot(time(win),ft(win),'b');mygridon;xlabel('time [s]')ylabel('force [N]')%axis([0, time(length(win)), min(f)-1, max(f)+1])%axis([0, win(length(win)), min(f)-1, 10])hold offmyreplot;subplot(2,2,3)plot(time(win),zt(win),'m');mygridon;xlabel('time [s]')ylabel('bristle displ. [m]')axis([0, time(max(size(win))), 0, max(abs(zt))])%axis([0, win(length(win)), 0, max(damping)])hold offmyreplot;subplot(2,2,4)plot(time(win),vrt(win),'r');mygridon;hold onplot(time(win),vbt(win),'b');xlabel('time [s]')ylabel('velocities [m/s]')axis([0, time(max(size(win))), 0, max(vrt)])hold offmyreplot;figure(2)title('Checking the nlfriction function')subplot(2,1,1)plot(time(win),countt(win));mygridon;xlabel('time [s]')ylabel('n.of iterations')%axis([0, time(length(win)), 0, max(...)])hold offmyreplot;subplot(2,1,2)plot(time(win),errt(win));mygridon;xlabel('time [s]')ylabel('abs. error')%axis([0, time(length(win)), 0, max(...)])hold off