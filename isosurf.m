function  isosurf(X,Y,Z,A)
p = patch(isosurface(X,Y,Z,A,0.7));
isonormals(X,Y,Z,A,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud