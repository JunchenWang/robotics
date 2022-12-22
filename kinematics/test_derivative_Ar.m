function test_derivative_Ar
r =  rand(3,1); 
dr = rand(3,1); 
dt = 1e-6; 
dA1 = (w_dr_A(r + dr * dt) - w_dr_A(r)) / dt; 
dA2 = derivative_Ar(r, dr); 
norm(dA1 - dA2)
