function [sv_list,vec0_list] = compute_deform_euclidean3d(x0,y0,z0,xf,yf,zf)
% x0,y0,xf,yf: (N1 x N2 x N3) double array
% x (y) - varies along rows (coloumns)
% Output 
% sv_list : N1 X N2 x N3 x 3: sv_list(:,1)< sv_list(:,2)
% vec0_list : N1 x N2 x N3 x 3 x 3: vec0_list(i,:,1) - vec0_list(i,:,2) -..
 
[N1,N2,N3] = size(x0); % Dimensions of the grid 

% Intialize the variables
sv_list = nan(N1,N2,N3,3); % Singular value of DF
vec0_list = nan(N1,N2,N3,3,3); % Eig. vector of DF at x0

for i = 2:N1-1
    for j = 2:N2-1
        for k = 2:N3-1
            dx = (x0(i+1,j,k)-x0(i-1,j,k))/2;
            dy = (y0(i,j+1,k)-y0(i,j-1,k))/2;
            dz = (z0(i,j,k+1)-z0(i,j,k-1))/2;

            DF11 = (xf(i+1,j,k)-xf(i-1,j,k))/(2*dx);
            DF12 = (xf(i,j+1,k)-xf(i,j-1,k))/(2*dy);
            DF13 = (xf(i,j,k+1)-xf(i,j,k-1))/(2*dz);

            DF21 = (yf(i+1,j,k)-yf(i-1,j,k))/(2*dx);
            DF22 = (yf(i,j+1,k)-yf(i,j-1,k))/(2*dy);
            DF23 = (yf(i,j,k+1)-yf(i,j,k-1))/(2*dz);

            DF31 = (zf(i+1,j,k)-zf(i-1,j,k))/(2*dx);
            DF32 = (zf(i,j+1,k)-zf(i,j-1,k))/(2*dy);
            DF33 = (zf(i,j,k+1)-zf(i,j,k-1))/(2*dz);
            
            DF = [DF11,DF12,DF13;DF21,DF22,DF23;DF31,DF32,DF33];
            [V,D] = eig(DF'*DF);
            [D,I] = sort(diag(D)); % Sort in ascending order  
        
            V0 = V(:,I); % Makes the matrix is V0: 2 x 2) :eigVec at x0
    
            sv_list(i,j,k,:) = sqrt(D); % Store singular values
            vec0_list(i,j,k,:,:) = V0;    % Store eigenvectors at x0
        end
    end 
end 

end 
