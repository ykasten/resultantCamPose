function [ v ] = hamiltonProduct( v1,v2 )
%HAMILTONPRODUCT calculate the Hamilton product of two quaternions
%v1 - the first quaternion
%v2 - the second quaternion
%v - the Hamilton product of the two quaternions
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019
a1=v1(1);
b1=v1(2);
c1=v1(3);
d1=v1(4);

a2=v2(1);
b2=v2(2);
c2=v2(3);
d2=v2(4);

v(1)=a1*a2-b1*b2-c1*c2-d1*d2;
v(2)=a1*b2+b1*a2+c1*d2-d1*c2;
v(3)=a1*c2-b1*d2+c1*a2+d1*b2;
v(4)=a1*d2+b1*c2-c1*b2+d1*a2;
end

