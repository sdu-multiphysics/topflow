function dphidg = dPHIdg(dx,dy,mu,alpha,dalpha,in6)
%dPHIdg
%    DPHIDG = dPHIdg(DX,DY,MU,ALPHA,DALPHA,IN6)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    26-May-2022 16:50:39

%Version: 0.9
u1 = in6(1,:);
u2 = in6(2,:);
u3 = in6(3,:);
u4 = in6(4,:);
u5 = in6(5,:);
u6 = in6(6,:);
u7 = in6(7,:);
u8 = in6(8,:);
dphidg = (dalpha.*dx.*dy.*(u1.*u3.*2.0+u1.*u5+u2.*u4.*2.0+u1.*u7.*2.0+u2.*u6+u3.*u5.*2.0+u2.*u8.*2.0+u3.*u7+u4.*u6.*2.0+u4.*u8+u5.*u7.*2.0+u6.*u8.*2.0+u1.^2.*2.0+u2.^2.*2.0+u3.^2.*2.0+u4.^2.*2.0+u5.^2.*2.0+u6.^2.*2.0+u7.^2.*2.0+u8.^2.*2.0))./3.6e+1;
