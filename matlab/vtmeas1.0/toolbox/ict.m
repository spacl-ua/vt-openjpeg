% ICT (RGB->YCbCr)
% Author: Han Oh, 3/7/2009
function out = ict(in)
a_R = 0.299;
a_G = 0.587;
a_B = 0.114;
out(:,:,1) = a_R*in(:,:,1)+a_G*in(:,:,2)+a_B*in(:,:,3);
out(:,:,2) = 0.5/(1-a_B)*(in(:,:,3)-out(:,:,1));
out(:,:,3) = 0.5/(1-a_R)*(in(:,:,1)-out(:,:,1));
end
