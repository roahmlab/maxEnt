mc = 10;  mp = 1;  l  = 0.5;  g = 9.81;

% Construct Equations
x0 = [ 0 0 0 0]';

x  = msspoly('x');
z  = msspoly('Tz');
[cbar,sbar] = unit_circle(z);
v  = msspoly('v',2);
u  = msspoly('u',1);


th0 = pi;
s   = -sbar;
c   = cbar;


H = [ mp*l^2 mp*l*c
      mp*l*c mc+mp  ];
C = [ 0 0 ; -mp*l*v(1)*s 0];
G = [ mp*g*l*s ; 0 ];
E = [ 0 ; 1];

F = E*u - C*v - G;

load cartpoleLQR;


detH = H(1,1)*H(2,2) - H(1,2)*H(2,1);
adjH = [ H(2,2) -H(1,2) ; -H(2,1) H(1,1) ];
Fcl  = subs(F,u,-K*[sbar;x;v]);

Sbar = S; Sbar(1,1) = 0;
V = S(1,1)*2*(1-cbar) + [sbar;x;v]'*Sbar*[sbar;x;v];



%[R,k] = pdecomp(V,z);


Vdot = diff(V,v)*adjH*Fcl + detH*diff(V,[z;x])*[ j*v(1)*z ; v(2)];

% Vdot = diff(V,[x;v])*Fcl(2:4) + ...
%        detH*R*((j*k*v(1)).*recomp(z,k,eye(length(k))));


ws = linspace(-10,10,100);
Ws = reshape(repmat(ws,200,1),[],1);
Ths = linspace(-pi,pi,200);
Zs = reshape(repmat(exp(j*Ths),1,100),[],1);

Vs = real(msubs(V,[z;x;v],[Zs 0*Zs Ws 0*Ws].'));

imagesc(Ths,Ws(1,:)',reshape(real(Vs),200,100)')
colorbar

% Randomized sampling test.

%% Now compute me a basin

prog = mssprog;

psd_to_herm = @(P,k) P(1:k,1:k) + P(k+(1:k),k+(1:k)) + j*(P(k+(1:k),1:k)-P(k+(1:k),1:k)');

Lphi = [ 1 ; z ; z^2; x ; v ; z*x ; z*v ];
[prog,GL] = new(prog,2*length(Lphi),'psd');
HL = psd_to_herm(GL,length(Lphi));
L = real(Lphi'*HL*Lphi);

rho = 0.5;

Vdphi = [mpmonomials([z;x;v],0:3);z^4];
[prog,GVd] = new(prog,2*length(Vdphi),'psd');
HVd = psd_to_herm(GVd,length(Vdphi));

[prog,gamma] = new(prog,1,'free');

sproc = -Vdot + L*(V-rho) + gamma*(x^2+v'*v+2*(1-cbar));
cond  = Vdphi'*HVd*Vdphi - sproc;
d = deg(cond,z);

[xx,pp,MM] = decomp(cond);
MM(:,any(full(pp)<0,2)) = [];
pp(any(full(pp)<0,2),:) = [];
cond = recomp(xx,pp,MM);

prog.eq = imag(cond);
prog.eq = real(cond);


% prog.eq = imag(z^d*cond);
% prog.eq = real(z^d*cond);

[prog,info] = sedumi(prog,gamma,0,struct());
prog(gamma)

       
surf(Ths,ws,reshape(real(Vs),200,100)'); shading flat
hold on;
contour(Ths,ws,reshape(real(Vs),200,100)',rho,'k*')
hold off


%% Pablo's method

prog = mssprog;

psd_to_herm = @(P,k) P(1:k,1:k) + P(k+(1:k),k+(1:k)) + j*(P(k+(1:k),1:k)-P(k+(1:k),1:k)');

Lphi = mpmonomials(z,0:2,x,0:2,v,0:1);
[prog,Lcoeff] = new(prog,length(Lphi),'free');
L = Lcoeff'*Lphi;

[prog,rho] = new(prog,1,'free');

Vdphi = [mpmonomials(z,0:5,x,0:4,v,0:4)];
[prog,GVd] = new(prog,2*length(Vdphi),'psd');
HVd = psd_to_herm(GVd,length(Vdphi));

[prog,gamma] = new(prog,1,'free');

sproc = L*Vdot + (s^2+v'*v+x^2)*(V-rho);
cond  = Vdphi'*HVd*Vdphi - sproc;
d = deg(cond,z);
prog.eq = imag(z^d*cond);
prog.eq = real(z^d*cond);

[prog,info] = sedumi(prog,-rho,0,struct());
prog(rho)
