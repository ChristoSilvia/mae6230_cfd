%mesh generation
% n = # cells, n+1 = # nodes
function C = alpha_mesh(alpha, y0, yn, n)
    
%find Delta y1
dy1 = (alpha-1)*(yn-y0)/(alpha^n-1);

%construct mesh
j = 1:(n-1);
C = dy1*(alpha.^j-1)./(alpha-1)+y0; 
C = [y0  C  yn];
end