function [vB] = Buu(B,u,v)

vB(:,1)=B(1)*u(:,2).*v(:,3);
vB(:,2)=B(2)*u(:,3).*v(:,1);
vB(:,3)=B(3)*u(:,1).*v(:,2);

end

