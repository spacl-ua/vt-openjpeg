% codes to start the VT tests
% Author: Feng Liu, the University of Arizona, 2/24/2017
set(0,'units','pixels')  
Pix_SS = get(0,'screensize');
d_space=min([floor(Pix_SS(3)-1536)/4;floor((Pix_SS(4)-512)/2);30]);
if(d_space<=0)
    d_space=10;
    sz=min([floor((Pix_SS(3)-d_space*4)/6)*2;floor((Pix_SS(3)-d_space*2)/2)*2]);
    warning('Your screen is too small. Instead of 512-by-512 stimuli, %d-by-%d stimuli will be used. This will affect your VT measurement results!',sz,sz);
end
cd([pwd '/' 'toolbox']);
comp=input('Component (0-Luminance, 1-Cb, 2-Cr): ');
cont=1;
while(cont==1)
    if(isempty(comp))
        comp=input('Component (0-Luminance, 1-Cb, 2-Cr): ');
    elseif((comp~=0)&&(comp~=1)&&(comp~=2))
        comp=input('Component (0-Luminance, 1-Cb, 2-Cr): ');
    else
        cont=0;
    end
end
if(comp==0)
    JPEG2000_3VT_lum; % Luminance Component
else
    JPEG2000_3VT_chrom; % Chrominance Component
end
cd ..;