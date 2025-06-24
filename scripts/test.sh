ls -1 "$1"/*region*gbk | while read line
do
	echo $line
	dir=$(echo $line | cut -d'/' -f1)
	#echo $dir
	file=$(echo $line | cut -d'/' -f3)
	echo $file
done
