
jnd=12800
echo "###########CMD options#########"
./demo_vts_kdu_format jnd=${jnd} extrap_method=3
echo ""
echo "##### REPLACE .cpp files###########"
./demo_vts_replace_cpp jnd=${jnd} extrap_method=3
