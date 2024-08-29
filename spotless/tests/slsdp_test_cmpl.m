prog = slsdp;
[prog,f1] = add(prog,10,'free');
[prog,f2] = add(prog,10,'free','complex');
[prog,l1] = add(prog,[5 3],'lor');
[prog,l2] = add(prog,[4 3],'lor','complex');