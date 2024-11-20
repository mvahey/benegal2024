function fit_dissociation(koff_fab,koff_IgG,B0)
%==========================================================================
% fit_dissociation(koff_fab,koff_IgG);
% Determines the crosslinking rate (kx) for the input kinetic parameters.
% Assumes that initial state consists of bivalently-bound IgG.
% koff_fab: dissociation rate for the Fab [1/s]
% koff_IgG: dissociation rate for the IgG [1/s]
% B0: Initial fraction of doubly-bound antibodies per HA (<=0.5)
%==========================================================================

% Fab binding parameters:
koff = koff_fab;     % Value for Fab
koff_eq = koff_IgG;  % Value for the IgG
dt = 0.01;

%==========================================================================
% First, generate the *effective* dissociation curves for IgG and Fab:
t = 0; A = B0; F = B0; T = 0;
a(1) = A; f(1) = F; t(1) = T;
n = 2;
while T<2000;
  A = a(n-1); F = f(n-1); T = t(n-1);
  dA = -koff_eq*A;
  dF = -koff*F;
  a(n) = A+dA*dt; f(n) = F+dF*dt; t(n) = T + dt;
  n = n + 1;
end
a_eq = a;
%==========================================================================


%==========================================================================
% Next, scan through kx values to find the best fit:
q = 1;
for kx = [[0:0.001:0.1] [0.11:0.01:1] [1.1:0.1:11]]
t = 0; A = 0; B = B0; T = 0;
a(1) = A; b(1) = B; t(1) = T;
n = 2;
while T<2000;
  A = a(n-1); B = b(n-1); T = t(n-1);
  dA = -koff*A-kx*A*(1-A-2*B)+koff*B;
  dB = kx*A*(1-A-2*B)-koff*B;
  a(n) = A+dA*dt; b(n) = B+dB*dt; t(n) = T + dt;
  n = n + 1;
end
%==========================================================================

plot(t,a_eq,'k',t,(a+b),'r'); 
ylim([0 B0]); drawnow; %pause;
E(q) = sum((a_eq-(a+b)).^2);
KX(q) = kx;
q = q + 1;
end

[val,ind] = min(E);

%==========================================================================
% Finally, plot the best fit value:
kx = KX(ind)
t = 0; A = 0; B = B0; T = 0;
a(1) = A; b(1) = B; t(1) = T;
n = 2;
while T<2000;
  A = a(n-1); B = b(n-1); T = t(n-1);
  dA = -koff*A-kx*A*(1-A-2*B)+koff*B;
  dB = kx*A*(1-A-2*B)-koff*B;
  a(n) = A+dA*dt; b(n) = B+dB*dt; t(n) = T + dt;
  n = n + 1;
end
plot(t,(a+b)/B0,'b'); hold on;
plot(t,a/B0,'r',t,b/B0,'r')
plot(t,a_eq/B0,'k'); xlim([0 2000]); ylim([0 1.1]); 
axis square;
drawnow;
hold off;