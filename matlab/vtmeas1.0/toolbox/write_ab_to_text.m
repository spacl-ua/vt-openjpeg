% Write VT measurement results to a text file
% Author: Yuzhang Lin, the University of Arizona, 3/5/2017
%Input
%       a: slopes,formated as:
%              [LH/HL1,HH1]
%              [LH/HL2,HH2]
%                  ...
%       b: intercept, same format as a
%       t: threshold for LL5
%       compno: component number. 1-Y; 2-Cb; 3-Cr.
%       VT index.
function write_ab_to_text(a,b,t,compno,idx_VT)



output_filename = [pwd,'/data/VT_table.txt'];

if ~exist(output_filename,'file')
   fid = fopen(output_filename,'w+');
   %write the column name at the first line
   fprintf(fid,'JND\tColor\tSubband\tThreshold\tStepsize\tSlope\n');
else
   %check if the firstline is the column name
   fid = fopen(output_filename,'a');
end
%% Write subbands threshold and slopes
% LL5 

str_band  = 'LL5';
threshold = t/256;
stepsize  = threshold;
slope     = 0;

fprintf(fid,'%d\t%d\t%s\t%.6f\t%.6f\t%.6f\n',...
              idx_VT,compno,str_band,threshold,stepsize,slope);

% HL5-->HH1
subband_name = {'LH','HL','HH'}

for level=5:-1:1
        
        for orient=1:3
                str_band  = [subband_name{orient},num2str(level)];
                
                if orient == 1 || orient==2 %HL/LH
                        threshold = b(level,1)/256; % normalized to [0,0.5]
                        slope     = a(level,1);
                        
                else %HH
                        threshold = b(level,2)/256; % normalized to [0,0.5]
                        slope     = a(level,2);
                end
                
                if (compno==1)
                        stepsize  = 1.03*threshold;
                else
                        stepsize  = threshold;
                end
                fprintf(fid,'%d\t%d\t%s\t%.6f\t%.6f\t%.6f\n',...
                        idx_VT,compno,str_band,threshold,stepsize,slope);
                
        end
end