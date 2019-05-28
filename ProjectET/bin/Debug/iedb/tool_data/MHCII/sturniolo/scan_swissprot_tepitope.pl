open(IN, "HLA_list.txt") or die "dd";
while($line=<IN>){
    %mat=();
    open(MAT, "$line") or die "$line";
    $ll=<MAT>;
    while($ll=<MAT>){
	@arr=split /\s+/, $ll;
	for($i=1; $i<=9; $i++){
	    $mat{$arr[0]}{$i}=$arr[$i];
	}
    }

    $oo="swissprot_".$line;
    open(OUT, ">$oo") or die "$oo";
    open(SW, "/home/pwang/db/uniprot_sprot.fasta") or die "dd";
    $counter=0;
    $fas="";
    while($ss=<SW>){
	if($ss=~/>/){
	    $end=length($fas)-9;
	    for($j=0; $j<=$end; $j++){
		$counter++;
		$score=0;
		$str=substr($fas, $j, 9);
		@str_arr=split //, $str;
		for($k=0; $k<=$#str_arr; $k++){
		    $index=$k+1;
		    $score=$score + $mat{$str_arr[$k]}{$index};
		}
		if($score<=-12.5){
		    $score=-12.5;
		}
		if($counter < 10000000){
		    print OUT "$str\t$score\n";
		}
	    }
	    $fas="";
	}
	else{
	    $ss=~/(\S+)/;
	    $fas=$fas.$1;
	}
    }
}
