function lwr = lwr_create(joints)
%Creates a KUKA LWR IV for Robotics Toolbox with DH parameters (without offset) 
%Returns robot struct with dynamic model for 7-DOF


    I_15= [0.003 0.003 0.003];
    I_6 = [0.002 0.002 0.002];
    I_7 = [0.001 0.001 0.001];
    Jm = 1e-5; %motor inertia
    G=[160 160 160 60 100 100 100]; %gear ratio
    B = 1e-5; %motor viscuous friction
    Tc = 1*[0.01 -0.01]; %link coulomb friction
    
    L1 = Link('d', 0.3105, 'a', 0, 'alpha', pi/2, 'm', 2, 'r', [0 0.1 0], 'I', I_15, 'Jm',Jm, 'G',G(1), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-170 170]));
    L2 = Link('d', 0, 'a', 0, 'alpha', -pi/2,     'm', 2, 'r', [0 0 0.1], 'I', I_15, 'Jm',Jm, 'G',G(2), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-120 120]));
    L3 = Link('d', 0.4, 'a', 0, 'alpha',-pi/2,    'm', 2, 'r', [0 0.1 0], 'I', I_15, 'Jm',Jm, 'G',G(3), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-170 170]));
    L4 = Link('d', 0, 'a', 0, 'alpha', pi/2,      'm', 2, 'r', [0 0 0.1], 'I', I_15, 'Jm',Jm, 'G',G(4), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-120 120]));
    L5 = Link('d', 0.39, 'a', 0, 'alpha', pi/2,   'm', 2, 'r', [0 0.1 0], 'I', I_15, 'Jm',Jm, 'G',G(5), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-170 170]));
    L6 = Link('d', 0, 'a', 0, 'alpha', -pi/2,     'm', 1, 'r', [0 0 0.04], 'I', I_6, 'Jm',Jm, 'G',G(6), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-120 120]));
    L7 = Link('d', 0.078, 'a', 0, 'alpha', 0,     'm', 0.1,'r',[0 0 0.00], 'I', I_7, 'Jm',Jm, 'G',G(7), 'B',B ,'Tc',Tc, 'qlim', deg2rad([-170 170]));
    
    lwr = SerialLink([L1 L2 L3 L4 L5 L6 L7], 'name', 'LWR 7DOF');

end
