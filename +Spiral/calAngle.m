function angle=calAngle(x0,y0,x1,y1,x2,y2)
% calculate the angle between two vectors 
A=[x1-x0,y1-y0];
B=[x2-x1,y2-y1];
y=det([A;B]);
x=dot(A,B);
angle=atan2(y,x);