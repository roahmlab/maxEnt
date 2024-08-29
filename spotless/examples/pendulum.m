z = msspoly('Tz'); % cos(th)+jsin(th)
[c,s] = unit_circle(z);
w = msspoly('w'); % thdot

psd_to_herm = @(P,k) P(1:k,1:k) + P(k+(1:k),k+(1:k)) + j*(P(k+(1:k),1:k)-P(k+(1:k),1:k)');

% Find a Lyapunov function proving stability of a damped pendulum
%
%  th.. = -th. - sin(th)
%
%
%  Search for a Lyapunov function:
%
%  V(th,th.) = q1*(1-cos(th)) + q2*sin(th)th. + q3*th.^2

prog = mssprog;

Vphi = [ 1 ; z ; z*w ; w];
[prog,GV] = new(prog,2*length(Vphi),'psd');
HV = psd_to_herm(GV,length(Vphi));
V = Vphi'*HV*Vphi;

prog.eq = trace(HV) - 1;

[R,k] = pdecomp(V,z);
Vdot = diff(V,w)*(-w-s) + R*((j*k*w).*recomp(z,k,eye(length(k))));

% Lphi = [ 1 ; z ; w ];
% [prog,GL] = new(prog,2*length(Lphi),'psd');
% HL = psd_to_herm(GL,length(Lphi));
% L = real(Lphi'*HL*Lphi);

Vdphi = [ 1 ; z ; z^2 ; z*w ; z^2*w ; w ];
[prog,GVd] = new(prog,2*length(Vdphi),'psd');
HVd = psd_to_herm(GVd,length(Vdphi));


[prog,gamma] = new(prog,1,'free');

sprog = prog;

%%  Vdot < -eps*(s^2+th.^2)
prog = sprog;

sproc = -Vdot + gamma*(s^2+w^2);

prog.eq = real(Vdphi'*HVd*Vdphi) - real(sproc);
prog.eq = imag(sproc);

[prog,info] = sedumi(prog,gamma,0,struct());

clean(subs(prog(V),z,c+j*s),1e-10)


%% Plot it

Ws = reshape(repmat(linspace(-5,5,400),200,1),[],1);
Ths = linspace(-pi,pi,200);
Zs = reshape(repmat(exp(j*Ths),1,400),[],1);

Vs = real(msubs(prog(V),[w;z],[Ws Zs].'));

imagesc(Ths,Ws(1,:)',reshape(real(Vs),200,400)')
colorbar


















