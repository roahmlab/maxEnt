% Simplest possible program

prog = slsdp;
x  = msspoly('x');
mn = monomials(x,0:10);
[prog,coeff] = prog.add(length(mn),'free');
f = coeff'*mn;

K = 100;
tt = linspace(-1,1,K);
yy = double(tt < 0);

err = yy - msubs(f,x,tt);

[prog,q] = add(prog,K+1,'lor');
t = q(1);
prog = constraint(prog,q(1+(1:K)) - err');

% [prog,t] = add_nn(prog,1);
% [prog,p] = add_nn(prog,K);
% [prog,m] = add_nn(prog,K);

% prog = constraint(prog, p' - (t - err));
% prog = constraint(prog, m' - (t + err));

prog = solve(prog,t);

fopt = prog.eval(f);
plot(tt,msubs(fopt,x,tt),tt,yy)

%% Test Nuclear norm calculation.
n=4;
k=2;
A = randn(n,k);

prog = slsdp;

[prog,P] = add(prog,n+k,'psd');
prog = constraint(prog,P(1:n,n+(1:k)) - A);

prog = solve(prog,trace(P));

%% Test a Hermitian matrix version these ideas.



