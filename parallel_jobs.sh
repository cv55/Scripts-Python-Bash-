#input files
sync_dir=$1
ne_file=$2

#output directory
out_dir="results"

start=`date +%s`
counter=0;

for sync in $sync_dir/*.sync 
do

	counter=$((counter+1));
	ne=$(cat $ne_file | awk -v counter2="$counter" 'NR==counter2')

	n=`basename $sync`
	echo "File under process: "$n;


	##keep pops file
	pops=$sync_dir/$n".pops"
	

	##make folder to store split files
	split_dir=$sync_dir"/splitdir"
	mkdir -p $split_dir


	##split file
	split -l 30000 $sync $split_dir"/"


	#generate *.pops file for each split file
	for a in $split_dir/a*
	do
		a_name=`basename $a`
		suff=".pops";
		`cp $pops $split_dir"/"$a_name$suff`

	done


	
	
	
	##run in parallel
	parallel 'python2 CLEAR.py --sync {1} --N {2}  --out {1}.out' ::: $split_dir/a[a-z] ::: $ne


	echo "Process finished for: "$n;

	mv $split_dir/*.out $out_dir/

	cd $out_dir
	python clear_output.py a*.out;

	awk '{if($1=="2L") print $0}' *.tsv > $n.clear

	rm a*.out
	rm *.tsv
	cd ../
	
	rm -r $split_dir
	


end=`date +%s`
runtime=$((end-start))
echo "Runtime was: " $runtime

done
