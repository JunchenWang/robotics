classdef LinearTrajectory
    properties
        tnodes
        Tds
        pps
        n
    end
    methods
        function obj = LinearTrajectory(tnodes, Tds)
            obj.tnodes = tnodes;
            obj.Tds = Tds;
            obj.n = length(tnodes) - 1;
            obj.pps = cell(1, obj.n);
            for i = 1 : obj.n
                [~,~,~,~,pp] = trapveltraj([0, 1], 2, 'EndTime', tnodes(i + 1) - tnodes(i));
                obj.pps{i} = pp{1};
            end
        end
        function [Td, vel, acc] = desired_pose(obj, t)
            flag = 0;
            for i = 1 : obj.n
                if obj.tnodes(i) <= t && t <= obj.tnodes(i + 1)
                    t = t - obj.tnodes(i);
                    Ts = obj.Tds{i};
                    Te = obj.Tds{i + 1};
                    pp = obj.pps{i};
                    ppd = fnder(pp, 1);
                    ppdd = fnder(pp, 2);
                    s = ppval(pp, t);
                    sd = ppval(ppd, t);
                    sdd = ppval(ppdd, t);
                    ps = Ts(1:3,4);
                    pe = Te(1:3,4);
                    Rs = Ts(1:3,1:3);
                    Re = Te(1:3,1:3);
                    pd = ps + (pe - ps)*s;
                    xe = logR(Rs'*Re)';
                    Rd = Rs * exp_w(xe*s);
                    Td = [Rd, pd; 0 0 0 1];
                    vel = [Rs * xe; (pe - ps)] * sd;
                    acc = [Rs * xe; (pe - ps)] * sdd;
                    flag = 1;
                    break;
                end
            end
            if flag == 0
                Td = obj.Tds(:,:,end);
                vel = zeros(6,1);
                acc = zeros(6,1);
            end
        end
    end
end