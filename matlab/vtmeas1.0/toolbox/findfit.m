% codes to linearly fit the VT versus the subband coefficient variance for
% luminance component
% Author: Feng Liu, the University of Arizona, 3/2/2017
function findfit(VTidx)
a=zeros(2,5);
b=zeros(2,5);
rsq=zeros(2,5);
loadhead={'LL','HL','LH','HH'};
for level=1:5
    %% read the variances and measured VTs for HL/LH
    
    sig21=zeros(4,1);
    vt1=zeros(4,1);
    for ki=1:4
        if(exist([pwd '/lumthresh/' loadhead{2} sprintf('%dvt%d_3afc_no_%d.mat',level,VTidx,ki)])==0)
            break;
        end
    end
    for kj=1:4
        if(exist([pwd '/lumthresh/' loadhead{3} sprintf('%dvt%d_3afc_no_%d.mat',level,VTidx,kj)])==0)
            break;
        end
    end
    if((ki~=4)||(kj~=4))
        cd ..;
        error('%d VT(s) for HL%d and %d VT(s) for LH%d VT(s) of luminance component at JND=%d are missing, take these measurements first.',5-ki,level,5-kj,level,VTidx);
    end
    for ki=1:4
        load([pwd '/lumthresh/' loadhead{2} sprintf('%dvt%d_3afc_no_%d.mat',level,VTidx,ki)]);
        sig20=sig2;
        t0=t;
        load([pwd '/lumthresh/' loadhead{3} sprintf('%dvt%d_3afc_no_%d.mat',level,VTidx,ki)]);
        if(t0<t)
            sig21(ki)=sig20;
            vt1(ki)=t0;
        else
            sig21(ki)=sig2;
            vt1(ki)=t;
        end
    end
    
    %% linearly fit for HL/LH
    
    [sig21,sortidx]=sort(sig21,'ascend');
    vt1=vt1(sortidx);
    p1=polyfit(sig21,vt1,1);
    yfit1=polyval(p1,sig21);
    yresid=vt1-yfit1;
    SSresid = sum(yresid.^2);
    SStotal = (length(vt1)-1) * var(vt1);
    rsq(2,level) = 1 - SSresid/SStotal;
    a(2,level)=p1(1);
    b(2,level)=p1(2);
    if(rsq(2,level)<0.9) % check if the measured VTs are linear enough v.s. the variances.
        warning('VT vs. coefficient variance for luminance HL/LH%d at JND=%d is off the linearity, remeasurement is highly recommended',level,VTidx);
    end
    
    %% read the variances and measured VTs for HH
    
    sig21=zeros(4,1);
    vt1=zeros(4,1);
    for ki=1:4
        if(exist([pwd '/lumthresh/' loadhead{4} sprintf('%dvt%d_3afc_no_%d.mat',level,VTidx,ki)])==0)
            cd ..;
            error('%d VT(s) for HH%d of luminance component at JND=%d are missing, take these measurements first.',5-ki,level,VTidx);
        end
        load([pwd '/lumthresh/' loadhead{4} sprintf('%dvt%d_3afc_no_%d.mat',level,VTidx,ki)]);
        sig21(ki)=sig2;
        vt1(ki)=t;
    end
    
    %% linearly fit for HH
    
    [sig21,sortidx]=sort(sig21,'ascend');
    vt1=vt1(sortidx);
    p1=polyfit(sig21,vt1,1);
    yfit1=polyval(p1,sig21);
    yresid=vt1-yfit1;
    SSresid = sum(yresid.^2);
    SStotal = (length(vt1)-1) * var(vt1);
    rsq(1,level) = 1 - SSresid/SStotal;
    a(1,level)=p1(1);
    b(1,level)=p1(2);
    if(rsq(1,level)<0.85) % check if the measured VTs are linear enough v.s. the variances.
        warning('VT vs. coefficient variance for luminance HH%d at JND=%d is off the linearity, remeasurement is highly recommended',level,VTidx);
    end
end
if(exist([pwd '/lumthresh/' sprintf('LL5vt%d_3afc_no_1.mat',VTidx)])~=0)
    load([pwd '/lumthresh/' sprintf('LL5vt%d_3afc_no_1.mat',VTidx)]);
    save([pwd '/data/' sprintf('LL5vt%d_3afc_no_1.mat',VTidx)],'t','sd','sig2');
else
    cd ..;
    error('1 VT(s) for LL5 of luminance component at JND=%d are missing, take these measurements first.',VTidx);
end
save([pwd '/data/' sprintf('abR_VT%d.mat',VTidx)],'a','b','rsq');
%write a,b to a text file for opj_compress encoder
write_ab_to_text(a,b,t,1,VTidx);