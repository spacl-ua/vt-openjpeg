% For Cr component, predict the future VTs based on those that have
% been measured by finding the linear coefficinets, assuming that the VT in
% a give subband is linear to the JND value.
% VTidx is the largest possible JND value that has been measured.
% Author: Feng Liu, the University of Arizona, 3/2/2017

function extrafit_Cr(VTidx)
loadhead={'LL','HL','LH','HH'};
bvec=zeros(2,5,VTidx);
VTidx_top=VTidx;
for ki=1:VTidx
    if(exist([pwd '/data/' sprintf('CrVT%d.mat',ki)])~=0)
        load([pwd '/data/' sprintf('CrVT%d.mat',ki)]);
        bvec(:,:,ki)=T(:,:);
    else
        VTidx_top=ki-1;
        break;
    end
end
ab=zeros(2,5);
bb=zeros(2,5);
rsq=zeros(2,5);
if(VTidx_top>1)
    bvec=bvec(:,:,1:VTidx_top);
    for ki=1:2
        for kj=1:5
            p1=polyfit((1:VTidx_top)',squeeze(bvec(ki,kj,:)),1);
            yfit1=polyval(p1,1:VTidx_top);
            yresid=squeeze(bvec(ki,kj,:))'-yfit1;
            SSresid = sum(yresid.^2);
            SStotal = (length(squeeze(bvec(ki,kj,:)))-1) * var(squeeze(bvec(ki,kj,:)));
            rsq(ki,kj) = 1 - SSresid/SStotal;
            ab(ki,kj)=p1(1);
            bb(ki,kj)=p1(2);
        end
    end
    save([pwd '/data/' sprintf('Cr_non_LL_abcoeff_up_to_%d.mat',VTidx_top)],'ab','bb');
end

bvec=zeros(1,VTidx);
VTidx_top=VTidx;
for ki=1:VTidx
    if(exist([pwd '/data/' loadhead{1} sprintf('5Crvt%d_3afc.mat',ki)])~=0);
        load([pwd '/data/' loadhead{1} sprintf('5Crvt%d_3afc.mat',ki)]);
        bvec(ki)=t;
    else
        VTidx_top=ki-1;
        break;
    end
end
if(VTidx_top>1)
    bvec=bvec(1:VTidx_top);
    p1=polyfit(1:VTidx_top,bvec,1);
    yfit1=polyval(p1,1:VTidx_top);
    yresid=bvec-yfit1;
    SSresid = sum(yresid.^2);
    SStotal = (length(bvec(:))-1) * var(bvec);
    rsq = 1 - SSresid/SStotal;
    ab=p1(1);
    bb=p1(2);
    save([pwd '/data/' sprintf('Cr_LL_abcoeff_up_to_%d.mat',VTidx_top)],'ab','bb');
end