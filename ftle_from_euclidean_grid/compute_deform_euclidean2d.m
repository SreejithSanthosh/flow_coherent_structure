function [sv_list,vec0_list] = compute_deform_euclidean2d(x0,y0,xf,yf)
% x0,y0,xf,yf: (N1 x N2) double array
% x (y) - varies along rows (coloumns)
% Output 
% sv_list : N1 X N2 x 2: sv_list(:,1)< sv_list(:,2)
% vec0_list : N1 x N2 x 2 x 2: vec0_list(i,:,1) - vec0_list(i,:,2) -..
 
[N1,N2] = size(x0); % Dimensions of the grid 

% Intialize the variables
sv_list = nan(N1,N2,2); % Singular value of DF
vec0_list = nan(N1,N2,2,2); % Eig. vector of DF at x0

for i = 2:N1-1
    for j = 2:N2-1
        
        DF11 = (xf(i,j+1)-xf(i,j-1))/(x0(i,j+1)-x0(i,j-1));
        DF12 = (xf(i+1,j)-xf(i-1,j))/(y0(i+1,j)-y0(i-1,j));
        DF21 = (yf(i,j+1)-yf(i,j-1))/(x0(i,j+1)-x0(i,j-1));
        DF22 = (yf(i+1,j)-yf(i-1,j))/(y0(i+1,j)-y0(i-1,j));
        
        DF = [DF11,DF12;DF21,DF22];

        [V,D] = eig(DF'*DF);
        [D,I] = sort(diag(D)); % Sort in ascending order  
    
        V0 = V(:,I); % Makes the matrix is V0: 2 x 2) :eigVec at x0

        sv_list(i,j,:) = sqrt(D); % Store singular values
        vec0_list(i,j,:,:) = V0;    % Store eigenvectors at x0
    end 
end 

end 
