while getopts "a:o:h" opt; do
	case $opt in
		a) alignfile=$OPTARG;;
		o) outp=$OPTARG;;
		h) printf "$helpmessage"
		   exit;;
		?) printf "$helpmessage"
                   exit;;
	esac
done

if [ $OPTIND -eq 1 ]; then
	printf "$helpmessage"
	exit
fi

outprefix=${outp}
outname=$outprefix"_vf_k"$kfilt"_I"$I"_X"$L"_filt-flag_filt-coord_scores.bed"
splitpref=$outprefix"_vf_k"$kfilt"_I"$I"_X"$L"_filt-flag_filt-coord_scores"


awk '$2==99 || $2==163 || $2==355 || $2==419{as = ""; ys = ""; for (i=12; i<=NF; i++){if ($i ~/^AS/){as = $i} if ($i ~/^YS/){ys = $i; break}} print $3"\t"$4"\t"$4+$9-1"\t"$1"\t"as"\t"ys}' $alignfile | awk '$2>0 && $3>$2 && $5~/^AS:i:/ && $6~/^YS:i:/{print}' > $outname

mkdir -p splits/

awk -v pref=$splitpref '{if ($4==storedval){a++; b[a] = $0; next} for (i in b){print b[i] > "splits/"pref"_map-"a".bed"} delete b; a = 1; b[a] = $0;; storedval = $4} END {for (i in b){print b[i] > "splits/"pref"_map-"a".bed"}}' $outname && gzip $outname
