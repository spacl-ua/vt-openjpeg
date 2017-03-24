
clear
%close all
suffix=''

folder_opj = '/home/sam/openjpeg-release/test_newest/jp2s_vt1'
folder_kdu = '/home/sam/visually_lossless_JPEG2000/ModifiedKakadu/bin/Linux-x86-64-gcc/jp2s_vt1'
list_opj=dir([folder_opj,'/*.bmp.jp2'])
list_kdu=dir([folder_kdu,'/*.bmp.jp2'])

num_files=numel(list_opj);

rate_vec_opj=zeros(1,num_files);
for i=1:num_files
    img_opj=imread([folder_kdu,'/',list_opj(i).name]);
    img_kdu=imread([folder_opj,'/',list_opj(i).name]);
    rate_vec_opj(i)=list_opj(i).bytes/numel(img_opj)*3;
    rate_vec_kdu(i)=list_kdu(i).bytes/numel(img_opj)*3;
    mse_vec(i)=sqrt(mean(double(img_opj(:))-double(img_kdu(:))));
end

fprintf('Average Bit Rate = %f (bpp) \n',mean(rate_vec_opj));

[~,idx_sort]=sort(rate_vec_kdu,'ascend');



figure
plot(rate_vec_kdu(idx_sort),mse_vec(idx_sort));
xlabel('rate');
ylabel('\Delta_{mse}')
title('difference in MSE')
figure
plot(rate_vec_kdu(idx_sort),abs(rate_vec_opj(idx_sort)-rate_vec_kdu(idx_sort)));
xlabel('rate');
ylabel('\Delta_{rate}')
title('difference in bit rate')
