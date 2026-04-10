io(1) = linio('IP/Constant', 1, 'input');
io(2) = linio('IP/inverted pendulum', 1, 'output');
op = operpoint('IP');
op.States.x = [0,0,0,0]';
option = linearizeOptions("StoreOffsets", "struct");
sys = linearize("IP", io, op, option);
