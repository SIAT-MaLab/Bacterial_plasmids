#usage: sh run.sh $fa $res_folder 
# $fa: sequences needed to be determined circular or not 
# $res_folder: the folder(uncreated or empty) for putting resluts
mkdir $2 
# For each input sequence, getting its begining 200bp and endding 200bp sequence. 
awk -v out_dir="$2" '{if(NR%2==1){a=$1}
else{print a"_1" > out_dir"/"a"_1.fa"; print substr($0,1,200) >> out_dir"/"a"_1.fa";
print a"_2"> out_dir"/"a"_2.fa";print substr($0,length-199) >> out_dir"/"a"_2.fa"}
}' $1 
# Entering working folder 
cd $2 
# Preparing parameters for running merger
rename 's/^>//' * 
ls *fa > dir 
sort dir | xargs -n 2 > dir2 
cut -d ' ' -f 1 dir2 > tmp1 
sed -i 's/$/_merger/' tmp1 
paste dir2 tmp1 > dir 
rm dir2 tmp1 
echo 'merger $1 $2 $3 /dev/null' > run_tmp 
# Running merger
cat dir | xargs -n 3 sh run_tmp 
rm run_tmp 
# Analyzing the results of merger. 
grep 'Identity' *merger > res1
sed -i 's/_1.fa_merger:# Identity://' res1
awk '{print $1,$2}' res1 > tmp1
awk -F '[/ ]' '{print $1,$2,$3,400-$3-$2}' tmp1 > tmp2
mv tmp2 res1
sort -k1,1 res1 > tmp1 
mv tmp1 res1
# res1 included 4 columns: 1. contig; 2. overalpping length; 3. length of the sequences(400bp) removing overlapping region; 4. mismatch number 
