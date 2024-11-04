%  An attempt to draw Photonic Band Gap meshes
%  for figures, but it's not that great.
%
% Last edit by Colin Joye, 7/2/03W

r = 0.1; % Circle radius: r/a < 0.5 for no overlap
a = 0.5; % Circle center separation
p = 20;  % number of points for the cylinder
h = 1;   % Cylinder height

G = [1 1 1 1 1 1 1 ;  % occupation matrix
     1 1 1 1 1 1 1 ;
     1 1 1 0 1 1 1 ;
     1 1 0 0 0 1 1 ;
     1 1 1 0 1 1 1 ;
     1 1 1 1 1 1 1 ;
     1 1 1 1 1 1 1 ];

figure(1)
clf(1)
opengl neverselect;
Lm = size(G,1);
Ln = size(G,2);
hold on;
for N=1:Ln,
    for M=1:Lm,
       [X,Y,Z] = cylinder([r r],p);
       pm = ((Lm+1)/2-M)*a;
       pn = ((Ln+1)/2-N)*a; 
       if( G(M,N)==1 ),
           surf( X+pm , Y+pn , Z*h, X );
       else
           line( X(1,:)+pm , Y(1,:)+pn , Z(1,:) );
       end
       colormap copper;
       shading interp; 
   end
end
grid on;
view([-27,76]);
hold off; 
axis equal;