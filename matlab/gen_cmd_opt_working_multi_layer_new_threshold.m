    % latest update: 
%       (1) both stepsize and thresholds are normlizaed to the range of
%       [0,0.5]
%       (2) use ":" as separater of two componets/layers
%% stepsize 
jnd_level=1;
para_file='cmd_option_multi_layer'

[a_lum,b_lum,b_c1,b_c2]=calc_VT_new_threshold(jnd_level);

b_lum(b_lum<128)=b_lum(b_lum<128);
b_c1(b_c1<128)=b_c1(b_c1<128);
b_c2(b_c2<128)=b_c2(b_c2<128);
%output from Feng's script is
%[LL5,HL5, LH5, HH5, HL4, LH4, HH4, HL3, LH3, HH3, HL2, LH2, HH2, HL1, LH1,HH1].
%Change the order to openjpeg friedly order:

str_stepsize=['-stepsize ',sprintf('%.4f,',[b_lum(1),b_lum(2:end-1)*1.03]/256),sprintf('%.4f:',b_lum(end)*1.03/256),...
    sprintf('%.4f,',[b_c1(1),b_c1(2:end-1)]/256),sprintf('%.4f:',b_c1(end)/256),...
    sprintf('%.4f,',[b_c2(1),b_c2(2:end-1)]/256),sprintf('%.4f',b_c2(end)/256)];

%% linear model slope and intercept 
str_slope=['-lin_model_slope '];
str_intercept=['-threshold '];
%[ HL1, LH1,HH1, ..., HL5, LH5,HH5,LL5];
%idx_reorder=[14,15,16,11,12,13,8,9,10,5,6,7,2,3,4,1]; %old order
%[LL5,HL5, LH5, HH5, HL4, LH4, HH4, HL3, LH3, HH3, HL2, LH2, HH2, HL1, LH1,HH1].
idx_reorder = 1:16;

for jnd_level=4:-1:1
    [a_lum,b_lum,b_c1,b_c2]=calc_VT_new_threshold(jnd_level);
    a_lum=a_lum(idx_reorder);
    b_lum=b_lum(idx_reorder)/256;
    b_c1=b_c1(idx_reorder)/256;
    b_c2=b_c2(idx_reorder)/256; %range 0-0.5
    
    
    str_slope=[str_slope,sprintf('%.4f,',a_lum(1:end-1)),sprintf('%.4f',a_lum(end))]
    str_intercept=[str_intercept,sprintf('%.4f,',b_lum(1:end-1)),sprintf('%.4f:',b_lum(end)),...
        sprintf('%.4f,',b_c1(1:end-1)),sprintf('%.4f:',b_c1(end)),...
        sprintf('%.4f,',b_c2(1:end-1)),sprintf('%.4f',b_c2(end))]
    
    %delimiter for VTs
    if jnd_level~=1
        str_slope=[str_slope,sprintf(':')];
        str_intercept=[str_intercept,sprintf(':')];
    end
end

fid=fopen(para_file,'w+');
fprintf(fid,'%s %s %s ',str_stepsize,str_slope,str_intercept);
fclose(fid)