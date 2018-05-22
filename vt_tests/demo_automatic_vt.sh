#!/bin/bash          
input_file=ki.bmp
#text file containing VTs
vt_file='../data/monitorA_VT1.txt'
opt_arg_stepsize="0.0027,0.0015,0.0015,0.0017,0.0018,0.0018,0.0020,0.0024,0.0024,0.0028,0.0035,0.0035,0.0071,0.0150,0.0150,0.0381:0.0028,0.0024,0.0024,0.0024,0.0039,0.0039,0.0148,0.0115,0.0115,0.0700,0.5000,0.5000,0.5000,0.5000,0.5000,0.5000:0.0026,0.0024,0.0024,0.0030,0.0037,0.0037,0.0046,0.0047,0.0047,0.0063,0.0082,0.0082,0.5000,0.5000,0.5000,0.5000"
opt_arg_slope="0.0000,0.0006,0.0006,0.0005,0.0004,0.0004,0.0002,0.0006,0.0006,0.0011,0.0015,0.0015,0.0020,0.0043,0.0043,0.0393"
opt_arg_threshold="0.0027,0.0015,0.0015,0.0016,0.0017,0.0017,0.0020,0.0024,0.0024,0.0027,0.0034,0.0034,0.0069,0.0145,0.0145,0.0370:0.0028,0.0024,0.0024,0.0024,0.0039,0.0039,0.0148,0.0115,0.0115,0.0700,0.5000,0.5000,0.5000,0.5000,0.5000,0.5000:0.0026,0.0024,0.0024,0.0030,0.0037,0.0037,0.0046,0.0047,0.0047,0.0063,0.0082,0.0082,0.5000,0.5000,0.5000,0.5000"

for jnd in {12800..12800}
do
#Using text file 
#../bin/opj_compress -i "${input_file}" -o "test_output_jnd${jnd}.jp2" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} 
#Decoding
#../bin/opj_decompress  -i "test_output_jnd${jnd}.jp2"  -o test_output_jnd${jnd}.jp2.bmp

#####################(FIXED STEPSIZE) use xx4s to replace xx5s ################
#../bin/opj_compress -i "${input_file}" -o "test_output_jnd${jnd}_extrap1.jp2" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} -extrap 1 
#../bin/opj_decompress  -i "test_output_jnd${jnd}_extrap1.jp2"  -o test_output_jnd${jnd}.jp2.bmp

#####################(FIXED STEPSIZE) probability summation pooling ################
#../bin/opj_compress -i "${input_file}" -o "test_output_jnd${jnd}_extrap2.jp2" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} -extrap 2 
#../bin/opj_decompress  -i "test_output_jnd${jnd}_extrap2.jp2"  -o test_output_jnd${jnd}.jp2.bmp

#####################(FIXED STEPSIZE) probability summation pooling with new interface -extrap AxBxC ################
#../bin/opj_compress -i "${input_file}" -o "test_output_jnd${jnd}_extrap2x0x0.jp2" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} -extrap 2x0x0 
#../bin/opj_decompress  -i "test_output_jnd${jnd}_extrap2x0x0.jp2"  -o test_output_jnd${jnd}.jp2.bmp
#####################(FIXED STEPSIZE) probability summation pooling with new interface -extrap AxBxC ################
jp2_file=${input_file}_jnd${jnd}_extrap2x0x0x3.jp2
mse_jp2_file=${input_file}_mse_jnd${jnd}_extrap2x0x0x3.jp2
#../bin/opj_compress -i "${input_file}" -o "${jp2_file}" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} -extrap 2x0x0x3
#../bin/opj_decompress  -i "${jp2_file}"  -o ${input_file}_jnd${jnd}_extrap2x0x0x3.jp2.bmp

#mse at jnd's bit rate
fsize=$(stat -c%s "${jp2_file}")
ratio=$(echo "scale=2;1024*1024*3/${fsize}" |bc)
echo "ratio=${ratio}"
../bin/opj_compress -i "${input_file}" -o "${mse_jp2_file}" -I -r ${ratio}
../bin/opj_decompress  -i "${mse_jp2_file}"  -o ${mse_jp2_file}.bmp

#####################(Dynamic STEPSIZE) probability summation pooling with new interface -extrap AxBxC ################
jp2_file=${input_file}_extrap2x0x0x3_dym.jp2;
mse_jp2_file=${input_file}_mse_extrap2x0x0x3_dym.jp2;
jp2_file_dec=${input_file}_extrap2x0x0x3_dym.jp2.bmp
../bin/opj_compress -i "${input_file}" -o "${jp2_file}" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} -extrap 2x0x0x3 -dynamic_stepsize 1
../bin/opj_decompress  -i "${jp2_file}"  -o ${jp2_file_dec}

fsize=$(stat -c%s "${jp2_file}")
ratio=$(echo "scale=2;1024*1024*3/${fsize}" |bc)
echo "ratio=${ratio}"
../bin/opj_compress -i "${input_file}" -o "${mse_jp2_file}" -I -r ${ratio}
../bin/opj_decompress  -i "${mse_jp2_file}"  -o ${mse_jp2_file}.bmp

#####################(DYNAMIC STEPSIZE)  ################
#../bin/opj_compress -i "${input_file}" -o "test_output_jnd${jnd}_dynamic_stepsize.jp2" -I -monitor ../data/asus-pa328q.txt -JND ${jnd} -dynamic_stepsize 1
#../bin/opj_decompress  -i "test_output_jnd${jnd}_dynamic_stepsize.jp2"  -o test_output_jnd${jnd}.jp2.bmp

done
