function [T] = transformation_matrix(l,x1,y1,x2,y2)

cost = (x2-x1)/l;
sint = (y2-y1)/l;

k = [cost sint 0;-sint cost 0;0 0 1];
y = zeros(3);

T = [k y;y k];

end