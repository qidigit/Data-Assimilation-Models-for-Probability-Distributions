function [vB] = Buu2(B,u,v)

vB1(1)=B(1)*u(2).*v(3);
vB1(2)=B(2)*u(3).*v(1);
vB1(3)=B(3)*u(1).*v(2);
vB=vB1';

end

