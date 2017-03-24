% Inverse ICT (YCbCr->RGB)
% Author: Han Oh, 3/7/2009
function out = invict(in)
out(:,:,1) = in(:,:,1) + 1.402*in(:,:,3);
out(:,:,2) = in(:,:,1) - 0.3441362862*in(:,:,2) - 0.7141362862*in(:,:,3);
out(:,:,3) = in(:,:,1) + 1.772*in(:,:,2);
end