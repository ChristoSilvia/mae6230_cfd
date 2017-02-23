%mesh_plot for part a

%build the mesh
C = alpha_mesh(1.05, 0, 1, 128);

%make it 2D for better visualization. 
%this is with equispaced points in x-direction
[xx, yy]=meshgrid(linspace(0,10,129), C);
mesh(xx,yy, xx.*0+1)
view(2)
colormap(bone)
export_fig -m2 -jpg -transparent 'meshgrid'

%extract the first 50 entries
[xx2, yy2] = meshgrid(linspace(0,10, 50), C(1:50)); 
mesh(xx2, yy2, xx2.*0+1)
view(2)
colormap(bone)
export_fig -m2 -jpg -transparent 'meshgrid_small'



