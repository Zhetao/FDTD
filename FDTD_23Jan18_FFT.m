clear
close
clc

%% **************************************************
% Define Simulation Parameters
dimen = 10000e-3; %Dimension of the whole system (waveguide)
X = 1e-3; %spatial sample period
dimD = round(dimen/X)+2; %Dimension [Nx] in nodes
c = zeros(1,dimD); %Sound speed in air
rho = zeros(1,dimD); %Air density
c0 = 343;
rho0 = 1.2;
cr = 0.5; %relative sound speed
c(1,:) = c0;
rho(1,:) = rho0;
K = 1./rho./((c*cr).^2);
Sc = c/c0*0.5; %Courant number

% T = X./c*Sc; %Time grid
T = X/c0*0.5; %Time grid
fs = 1/T; %Sample rate (Hz)
simTime = 30e-3; %simulation time, in seconds

alpha = 1e-3; %attenuation in PML
maxN = round(simTime/T);
maxc = 1.5; %caxis limits

freq = 1500;      %source frequency
% n0 = round(fs/freq/2);    % Initial delay (samples)
n=0:maxN;
srcFn=sin(2*pi*freq*(n-2)*T);
srcFn(1:2) = 0;
% srcFn(569:end) = 0;

%mod params
fmod = 400;  L = 1e-2; t0 = L/c0;
t1 = 0.77/2/freq/2;
xmod1 = 100;
xmod2 = xmod1+214+9;
recPos = xmod2+50;
%% **************************************************
% Allocate Matrices
p = zeros(2,dimD); %Pressure matrix (Nx Ny), (:,:,1) is n+1, (:,:,2) is n
ux = zeros(2,dimD); %velocity x matrix (Nx Ny)
att = zeros(1,dimD); %attenuation matrix
ll = 2:dimD-1;
pres = zeros(maxN,dimD);
ures = zeros(maxN,dimD);
xrec = zeros(1,maxN);
xrec = zeros(1,maxN);


att(2:20) = alpha;
att(end-19:end-1) = alpha;
%Visulize
% pl = p(1,:);
% h2 = plot(pl);
% ylim([-maxc maxc]);
% xlim([2 dimD-1]);
%% **************************************************
% Prepare Simulation
tic
for nn = 2:maxN-1
    %source here
    p(2,3) = srcFn(nn);
    
    %time-modulated phase
    c(xmod1:xmod1+9) = L/(t0+t1+t1*sin(2*pi*fmod*nn*T));
    rho(xmod1:xmod1+9) = (c0/(L/(t0+t1+t1*sin(2*pi*fmod*nn*T))))*rho0;
%     rhorec(nn) = rho(100);
%     crec(nn) = c(100);
    
    %time-modulated phase2
    c(xmod2:xmod2+9) = L/(t0+t1+t1*cos(2*pi*fmod*nn*T));
    rho(xmod2:xmod2+9) = (c0/(L/(t0+t1+t1*cos(2*pi*fmod*nn*T))))*rho0;
    
    %update grid
    ux(1,ll) = ux(2,ll)- T./rho(ll)/X.*(p(2,ll+1)-p(2,ll));
    %p(1,ll) = p(2,ll) -  T/X*c(2:dimD-1).^2.*(rho(ll).*ux(1,ll)-rho(ll-1).*ux(1,ll-1)); 
    p(1,ll) = p(2,ll) -  T/X*c(2:dimD-1).^2.*rho(ll).*(ux(1,ll)-ux(1,ll-1)); 

    ux(2,ll) = ux(1,ll);
    p(2,ll) = p(1,ll);
    
%     pres(nn,2:dimD-1) = p(2,ll);
%     ures(nn,2:dimD-1) = ux(2,ll);    
    xrec(1,nn) = p(2, recPos);
%     
%     %plot
     if mod(nn, 100) == 0
         display(nn);
%         pl = p(1,:);
%         plot(pl);
%         ylim([-maxc maxc]);
%         xlim([2 dimD]);
%     %     set(h2,'yData',pl);
%         title(sprintf('Sample number %d out of %d\n',nn,maxN));
%     %     pause(0.1)
%         drawnow;   
     end
end
toc
%
% figure(1)
% plot(xrec)
% 
% figure(2)
indexList = 1:maxN;
fList = indexList*fs/maxN;
plot(fList,abs(fft(xrec)))
xlim([0 3000])
% figure(1)
% plot(rhorec)
% figure(2)
% plot(crec)