% Collect the response given a stimulus at a certain strength
% Author: Feng Liu, the University of Arizona, 2/24/2017

function [response, err] = JND_trial_3VT(B, comp, level, band, blk_size_m, blk_size_n, delta_base, delta, sig2)

ofst=0;

set(0,'units','pixels')  
Pix_SS = get(0,'screensize');
d_space=min([floor(Pix_SS(3)-1536)/4;floor((Pix_SS(4)-512)/2);30]);
if(d_space<=0)
    d_space=10;
    sz=min([floor((Pix_SS(3)-d_space*4)/6)*2;floor((Pix_SS(3)-d_space*2)/2)*2]);
    startp=(512-sz)/2;
else
    sz=512;
    startp=0;
end
if comp ~= 0
    rgb_rx = 2^(B-1)*ones(sz,sz,3);   %pure gray color image
    rgb_x = 2^(B-1)*ones(sz,sz,3);   %pure gray color image
    ycbcr_rx = zeros(sz,sz,3);
    ycbcr_x = zeros(sz,sz,3);
end

mu = 0; %mean of codeblock

x = 2^(B-1)*ones(sz,sz);
y = zeros(sz,sz);

org_y = y;
idx = -3*level + band + 16;
[m, n] = size(y);
m = m*2^(-level);
n = n*2^(-level);
if band == 0
    y0 = y(1:m,1:m);
elseif band == 1
    y0 = y(1:m,m+1:2*m);
elseif band == 2
    y0 = y(m+1:2*m,1:m);
else
    y0 = y(m+1:2*m,m+1:2*m);
end

org_y0 = y0;


blk_size_m = min(blk_size_m, m);
blk_size_n = min(blk_size_n, n);
num_blk_m = m/blk_size_m;
num_blk_n = n/blk_size_n;

%% First let's generate the "background image"

uni_blk = -sqrt(12*sig2)/2+sqrt(12*sig2)*rand(blk_size_m,blk_size_n);
laplacian_blk = randl(blk_size_m,blk_size_n,mu,sqrt(sig2));

if(delta_base~=0)
    %noise shaping
    if (band == 0) %LL subband
        q_blk_base = uni_blk; %Uniform wavelet coefficients
        ii = find(q_blk_base>delta_base | q_blk_base<-delta_base ); %find non-deadzone coefficients
        coeffs=fix(q_blk_base(ii)/delta_base);
        q_blk_base(ii) = q_blk_base(ii)-(sign(coeffs).*(abs(coeffs)+0.5)*delta_base);
    else %HL/LH/HH subband
        q_blk_base = laplacian_blk; %generate wavelet coefficients according to Laplacian distribution
        ii = find(q_blk_base>delta_base | q_blk_base<-delta_base); %find non-deadzone coefficients
        coeffs=fix(q_blk_base(ii)/delta_base);
        q_blk_base(ii) = q_blk_base(ii)-(sign(coeffs).*(abs(coeffs)+0.5)*delta_base);
    end
    
    %center positioning
    y0(m/2-blk_size_m/2+1:m/2+blk_size_m/2,n/2-blk_size_n/2+1:n/2+blk_size_n/2) = org_y0(m/2-blk_size_m/2+1:m/2+blk_size_m/2,n/2-blk_size_n/2+1:n/2+blk_size_n/2)+q_blk_base;
    
    %inverse wavelet transform
    if band == 0
        y(1:m,1:m) = y0;
    elseif band == 1
        y(1:m,m+1:2*m) = y0;
    elseif band == 2
        y(m+1:2*m,1:m) = y0;
    else
        y(m+1:2*m,m+1:2*m) = y0;
    end
    
    if (level == 0)
        x = y;
    else
        try
            x = lwt(y,-level);
        catch
            mex lwt.cpp;
            x = lwt(y,-level);
        end
    end
else
    x=zeros(sz,sz);
end

%% Next, generate the image with additional stimulus
if(delta~=0)
    %noise shaping
    if (band == 0) %LL subband
        q_blk = uni_blk; %Uniform wavelet coefficients
        ii = find(q_blk>delta | q_blk<-delta ); %find non-deadzone coefficients
        coeffs=fix(q_blk(ii)/delta);
        q_blk(ii) = q_blk(ii)-(sign(coeffs).*(abs(coeffs)+0.5)*delta);
    else %HL/LH/HH subband
        q_blk = laplacian_blk; %generate wavelet coefficients according to Laplacian distribution
        ii = find(q_blk>delta | q_blk<-delta); %find non-deadzone coefficients
        coeffs=fix(q_blk(ii)/delta);
        q_blk(ii) = q_blk(ii)-(sign(coeffs).*(abs(coeffs)+0.5)*delta);
    end
    
    %center positioning
    y0(m/2-blk_size_m/2+1:m/2+blk_size_m/2,n/2-blk_size_n/2+1:n/2+blk_size_n/2) = org_y0(m/2-blk_size_m/2+1:m/2+blk_size_m/2,n/2-blk_size_n/2+1:n/2+blk_size_n/2)+q_blk;
    
    %inverse wavelet transform
    if band == 0
        y(1:m,1:m) = y0;
    elseif band == 1
        y(1:m,m+1:2*m) = y0;
    elseif band == 2
        y(m+1:2*m,1:m) = y0;
    else
        y(m+1:2*m,m+1:2*m) = y0;
    end
    
    if (level == 0)
        rx = y;
    else
        try
            rx = lwt(y,-level);
        catch
            mex lwt.cpp;
            rx = lwt(y,-level);
        end
    end
else
    rx=zeros(sz,sz);
end

%display x and rx randomly
stimuli_show = randi(3);
if comp ~= 0  %truecolor image
    big_img = uint8(zeros(sz+d_space*2+ofst,sz*3+d_space*4,3));
    ycbcr_rx(:,:,comp+1) = rx;
    rgb_rx = (invict(ycbcr_rx)+2^(B-1));
    
    ycbcr_x(:,:,comp+1) = x;
    rgb_x = (invict(ycbcr_x)+2^(B-1));
    
    if (stimuli_show == 1)
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space+1:d_space+sz,:) = uint8(rgb_rx(startp+1:startp+sz,startp+1:startp+sz,:));
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*2+sz+1:d_space*2+sz*2,:) = uint8(rgb_x(startp+1:startp+sz,startp+1:startp+sz,:));
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*3+sz*2+1:d_space*3+sz*3,:) = uint8(rgb_x(startp+1:startp+sz,startp+1:startp+sz,:));
    elseif (stimuli_show == 2)
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space+1:d_space+sz,:) = uint8(rgb_x(startp+1:startp+sz,startp+1:startp+sz,:));
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*2+sz+1:d_space*2+sz*2,:) = uint8(rgb_rx(startp+1:startp+sz,startp+1:startp+sz,:));
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*3+sz*2+1:d_space*3+sz*3,:) = uint8(rgb_x(startp+1:startp+sz,startp+1:startp+sz,:));
    elseif (stimuli_show == 3)
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space+1:d_space+sz,:) = uint8(rgb_x(startp+1:startp+sz,startp+1:startp+sz,:));
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*2+sz+1:d_space*2+sz*2,:) = uint8(rgb_x(startp+1:startp+sz,startp+1:startp+sz,:));
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*3+sz*2+1:d_space*3+sz*3,:) = uint8(rgb_rx(startp+1:startp+sz,startp+1:startp+sz,:));
    end
else
    big_img = zeros(sz+d_space*2+ofst,sz*3+d_space*4);
    rx = rx + 2^(B-1);
    x = x + 2^(B-1);
    if (stimuli_show == 1)
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space+1:d_space+sz) = rx(startp+1:startp+sz,startp+1:startp+sz);
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*2+sz+1:d_space*2+sz*2) = x(startp+1:startp+sz,startp+1:startp+sz);
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*3+sz*2+1:d_space*3+sz*3) = x(startp+1:startp+sz,startp+1:startp+sz);
    elseif (stimuli_show == 2)
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space+1:d_space+sz) = x(startp+1:startp+sz,startp+1:startp+sz);
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*2+sz+1:d_space*2+sz*2) = rx(startp+1:startp+sz,startp+1:startp+sz);
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*3+sz*2+1:d_space*3+sz*3) = x(startp+1:startp+sz,startp+1:startp+sz);
    elseif (stimuli_show == 3)
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space+1:d_space+sz) = x(startp+1:startp+sz,startp+1:startp+sz);
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*2+sz+1:d_space*2+sz*2) = x(startp+1:startp+sz,startp+1:startp+sz);
        big_img(d_space+1+ofst:d_space+sz+ofst,d_space*3+sz*2+1:d_space*3+sz*3) = rx(startp+1:startp+sz,startp+1:startp+sz);
    end
end


% set up the timer
tt = timer ;
tt.timerfcn = 'uiresume' ;
tt.startdelay = 10 ;  %displaying time (sec)

callstr = ['set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))); uiresume; '] ;
%display
if comp ~= 0
    fh = figure('keypressfcn',callstr);imshow(big_img,'InitialMagnification',100,'Border','tight');
else
    fh = figure('keypressfcn',callstr);imshow(big_img,'InitialMagnification',100,'DisplayRange',[0 2^B-1],'Border','tight');
end
start(tt) ;
uiwait ;
disp('Where is the signal (1,2,3)?');
reply = char(get(fh,'Userdata'));
stop(tt) ;
delete(tt) ;
close(fh) ;

cont = 1;
while (cont)
    if reply == '1'
        if stimuli_show == 1
            response = 1;  %correct
        else
            response = 0;  %incorrect
        end
        cont = 0;
    elseif reply == '2'
        if stimuli_show == 2
            response = 1;   %correct
        else
            response = 0;   %incorrect
        end
        cont = 0;
    elseif reply == '3'
        if stimuli_show == 3
            response = 1;   %correct
        else
            response =0;   %incorrect
        end
        cont=0;
        
    else
        reply = input('Where is the signal (1,2,3)? ','s');
    end
end


%Quantization distortion information
abs_dist_imgx = abs(rx-x);
abs_dist_imgy = abs(org_y-y);

err.xabs_max = max(abs_dist_imgx(:));
err.yabs_max = max(abs_dist_imgy(:));
err.delta = delta;
err.xabs_mean = mean(abs_dist_imgx(:));
err.yabs_mean = mean(abs_dist_imgx(:));
ind = find(abs_dist_imgy(:)>delta/2);
err.yabs_big_mean = mean(abs_dist_imgy(ind));



