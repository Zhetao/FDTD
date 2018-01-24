% This program describes a moving 1-D wave
% using the finite difference method
clc
close all;
clear all;


%-------------------------------------------------------------------------%
% Initialization

c = 343;
rho = 1.2;

Nx = 402;       % x-Grids
dx = 1e-3;         % Step size
x(:,1) = (0:Nx-1)*dx;

                
T = 2058;       % Total number of time steps
f = 1500;         % frequency of source

dt = dx/c*0.5;     % Time-Step
t(:,1)= (0:T-1)*dt;

% v = 500;        % Wave velocity
% c = v*(dt/dx);  % CFL condition

P = zeros(T,Nx);  % P(t,x) = P(time,space)
U = zeros(T,Nx);

s1 = T;  
%-------------------------------------------------------------------------%

% Initial condition

P((1:s1),3) = sin(2*pi*f.*t(1:s1));
%P((1:s1),4) = sin(2*pi*f.*t(1:s1));
%-------------------------------------------------------------------------%
%%
% Finite Difference Scheme
for j = 2:T
    P((1:s1),3) = sin(2*pi*f.*t(1:s1));
    for i = 2:Nx-1
        U(j,i) = U(j-1,i) - 1.0/rho*dt/dx*(P(j-1,i+1)-P(j-1,i));%why i+1 not i
        P(j,i) = P(j-1,i) - rho*c*c*dt/dx*(U(j,i)-U(j,i-1));
        if mod(j,20) == 0
            display(U(j,i))
        end
    end 
    
end
%-------------------------------------------------------------------------%
% Plot the results
%{
plot_times = [50 150 250 310];

for i = 1:4
    
  figure(i)
  k = plot_times(i);
  plot(x,U(k,:),'linewidth',2);
  grid on;
  axis([min(x) max(x) -2 2]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(t(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize') 
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end
%}
%-------------------------------------------------------------------------%
% Movie for the travelling wave

for j = 1:T              
  plot(x,P(j,:),'linewidth',2);
  grid on;
  axis([min(x) max(x) -2 2]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize') 
  set(h,'FontSize',14);
  fh = figure(5);
  set(fh, 'color', 'white'); 
  F=getframe;
            
end

movie(F,T,1)


%-------------------------------------------------------------------------%

