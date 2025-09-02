function [ codebook ] = Codebook_Generation

%% Generation Setup (downlink)
% -------------- Construct the generating vector -------------------------
u0 = [1 -1 -1 -1].';
u1 = [1 -1i 1 1i].';
u2 = [1 1 -1 1].';
u3 = [1 1i 1 -1i].';
u4 = [1 (-1-1i)/sqrt(2) -1i (1-1i)/sqrt(2)].';
u5 = [1 (1-1i)/sqrt(2) 1i (-1-1i)/sqrt(2)].';
u6 = [1 (1+1i)/sqrt(2) -1i (-1+1i)/sqrt(2)].';
u7 = [1 (-1+1i)/sqrt(2) 1i (1+1i)/sqrt(2)].';

u = [u0 u1 u2 u3 u4 u5 u6 u7];
% -------------- Construct the Codebook -------------------------
for n = 1:8
    W(:,:,n) = eye(4)-2*u(:,n)*u(:,n)'/(u(:,n)'*u(:,n));
end

codebook(:,:,1) = W(:,:,1)/2;
codebook(:,:,2) = W(:,:,2)/2;
codebook(:,:,3) = W(:,[3 2 1 4],3)/2;
codebook(:,:,4) = W(:,[3 2 1 4],4)/2;
codebook(:,:,5) = W(:,:,5)/2;
codebook(:,:,6) = W(:,:,6)/2;
codebook(:,:,7) = W(:,[1 3 2 4],7)/2;
codebook(:,:,8) = W(:,[1 3 2 4],8)/2;

end
