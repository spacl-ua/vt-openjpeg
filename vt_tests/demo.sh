#!/bin/bash          
input_file=test.bmp
#text file containing VTs
vt_file='../data/monitorA_VT1.txt'
opt_arg_stepsize="0.0027,0.0015,0.0015,0.0017,0.0018,0.0018,0.0020,0.0024,0.0024,0.0028,0.0035,0.0035,0.0071,0.0150,0.0150,0.0381:0.0028,0.0024,0.0024,0.0024,0.0039,0.0039,0.0148,0.0115,0.0115,0.0700,0.5000,0.5000,0.5000,0.5000,0.5000,0.5000:0.0026,0.0024,0.0024,0.0030,0.0037,0.0037,0.0046,0.0047,0.0047,0.0063,0.0082,0.0082,0.5000,0.5000,0.5000,0.5000"
opt_arg_slope="0.0000,0.0006,0.0006,0.0005,0.0004,0.0004,0.0002,0.0006,0.0006,0.0011,0.0015,0.0015,0.0020,0.0043,0.0043,0.0393"
opt_arg_threshold="0.0027,0.0015,0.0015,0.0016,0.0017,0.0017,0.0020,0.0024,0.0024,0.0027,0.0034,0.0034,0.0069,0.0145,0.0145,0.0370:0.0028,0.0024,0.0024,0.0024,0.0039,0.0039,0.0148,0.0115,0.0115,0.0700,0.5000,0.5000,0.5000,0.5000,0.5000,0.5000:0.0026,0.0024,0.0024,0.0030,0.0037,0.0037,0.0046,0.0047,0.0047,0.0063,0.0082,0.0082,0.5000,0.5000,0.5000,0.5000"
#Using text file 
../bin/opj_compress -i "${input_file}" -o "test_output.jp2" -I -vt ${vt_file}
../bin/opj_decompress -i test_output.jp2 -o test_output.jp2.bmp 
#Default stepsize, no linear model
../bin/opj_compress -i "${input_file}" -o "test_output_vt_only.jp2" -I -threshold ${opt_arg_threshold}
../bin/opj_decompress -i test_output_vt_only.jp2 -o test_output_vt_only.jp2.bmp 
#custom stepsize, with linear model
../bin/opj_compress -i "${input_file}" -o "test_output_custom_stepsize_linear_model.jp2" -I -stepsize ${opt_arg_stepsize} -lin_model_slope ${opt_arg_slope} -threshold ${opt_arg_threshold}
../bin/opj_decompress -i test_output_custom_stepsize_linear_model.jp2 -o test_output_custom_stepsize_linear_model.jp2.bmp 
../bin/opj_decompress -i test_output_custom_stepsize_linear_model.jp2 -o test_output.bmp 
#custom stepsize, without linear model
../bin/opj_compress -i "${input_file}" -o "test_output_custom_stepsize_no_linear_model.jp2" -I -stepsize ${opt_arg_stepsize}  -threshold ${opt_arg_threshold}
../bin/opj_decompress -i test_output_custom_stepsize_no_linear_model.jp2 -o test_output_custom_stepsize_no_linear_model.jp2.bmp 

#####################New functions for v2.0: compression at given JND ################


#####   Interface for releasing : most basic interface - Just JND   #####
target_jnd=1
jp2_file=${input_file}_jnd${target_jnd}.jp2
mse_jp2_file=${input_file}_mse_jnd${target_jnd}.jp2
../bin/opj_compress -i "${input_file}" -o "${jp2_file}" -I -monitor ../data/asus-pa328q.txt -JND ${target_jnd}
../bin/opj_decompress  -i "${jp2_file}"  -o ${jp2_file}.bmp


#####   Interface for releasing : finer tunning   #####

jp2_file=${input_file}_jnd${target_jnd}_tune.jp2
mse_jp2_file=${input_file}_mse_jnd${target_jnd}_tune.jp2
../bin/opj_compress -i "${input_file}" -o "${jp2_file}" -I -monitor ../data/asus-pa328q.txt -JND ${target_jnd} -perceptual_rate_alloc 2 -lum_chrom_balance 3
../bin/opj_decompress  -i "${jp2_file}"  -o ${jp2_file}.bmp

exit




#####   Interface for development  #####


jp2_file=${input_file}_jnd${target_jnd}_extrap2x0x0x3.jp2
mse_jp2_file=${input_file}_mse_jnd${target_jnd}_extrap2x0x0x3.jp2

../bin/opj_compress -i "${input_file}" -o "${jp2_file}" -I -monitor ../data/asus-pa328q.txt -JND ${target_jnd} -extrap 2x0x0x3
../bin/opj_decompress  -i "${jp2_file}"  -o ${input_file}_jnd${target_jnd}_extrap2x0x0x3.jp2.bmp
#mse at jnd's bit rate
fsize=$(stat -c%s "${jp2_file}")
ratio=$(echo "scale=2;768*512*3/${fsize}" |bc)
echo "ratio=${ratio}"
../bin/opj_compress -i "${input_file}" -o "${mse_jp2_file}" -I -r ${ratio}
../bin/opj_decompress  -i "${mse_jp2_file}"  -o ${mse_jp2_file}.bmp


