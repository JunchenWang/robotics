function R = RPY2R(rpy)
R = EulZYX2R([rpy(3), rpy(2), rpy(1)]);