function ret = Q(r)
ret = [r(1) -r(2) -r(3) -r(4);
     r(2)  r(1) -r(4)  r(3);
     r(3)  r(4)  r(1) -r(2);
     r(4) -r(3)  r(2)  r(1)];
    