% Put the non-LL subband measurement of Cb VTs at JND level of VTidx into the format for
% further synthesis.
% Author: Feng Liu, the University of Arizona, 3/2/2017

function findCbVT(VTidx)
loadhead={'LL','HL','LH','HH'};
T=zeros(2,5);
for ki=1:5
    if(exist([pwd '/chromthresh/' loadhead{2} sprintf('%dCbvt%d_3afc.mat',ki,VTidx)])~=0)
        load([pwd '/chromthresh/' loadhead{2} sprintf('%dCbvt%d_3afc.mat',ki,VTidx)]);
        t0=t;
    else
        cd ..;
        error(['Please measure the VT for Cb ' loadhead{2} sprintf('%d subband at JND=%d first.',ki,VTidx)]);
    end
    if(exist([pwd '/chromthresh/' loadhead{3} sprintf('%dCbvt%d_3afc.mat',ki,VTidx)])~=0)
        load([pwd '/chromthresh/' loadhead{3} sprintf('%dCbvt%d_3afc.mat',ki,VTidx)]);
        if(t>t0);
            t=t0;
        end
        T(2,ki)=t;
    else
        cd ..;
        error(['Please measure the VT for Cb ' loadhead{3} sprintf('%d subband at JND=%d first.',ki,VTidx)]);
    end
    if(exist([pwd '/chromthresh/' loadhead{4} sprintf('%dCbvt%d_3afc.mat',ki,VTidx)])~=0)
        load([pwd '/chromthresh/' loadhead{4} sprintf('%dCbvt%d_3afc.mat',ki,VTidx)]);
        T(1,ki)=t;
    else
        cd ..;
        error(['Please measure the VT for Cb ' loadhead{4} sprintf('%d subband at JND=%d first.',ki,VTidx)]);
    end
end
if(exist([pwd '/chromthresh/' loadhead{1} sprintf('5Cbvt%d_3afc.mat',VTidx)])~=0)
    load([pwd '/chromthresh/' loadhead{1} sprintf('5Cbvt%d_3afc.mat',VTidx)]);
    save([pwd '/data/' loadhead{1} sprintf('5Cbvt%d_3afc.mat',VTidx)],'t','sd');
else
    cd ..;
    error(['Please measure the VT for Cb ' sprintf('LL5 subband at JND=%d first.',VTidx)]);
end
save([pwd '/data/' sprintf('CbVT%d.mat',VTidx)],'T');
%write a,b to a text file for opj_compress encoder
write_ab_to_text(zeros(5,2),T,t,2,VTidx);