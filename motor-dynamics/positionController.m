function u = positionController(t, y, yd, K, velController)
vd = @(t) velocity_command(t, y, yd, K);
u = velController(t, y, vd);


    function vd = velocity_command(t, y, yd, K)
        ydt = yd(t);
        vd = [K * (ydt(1) - y(1)) + ydt(2); 0];
    end

end

