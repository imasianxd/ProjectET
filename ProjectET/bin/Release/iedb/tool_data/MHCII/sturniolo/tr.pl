open(IN, "list.txt") or die "dd";
while($line=<IN>){
    $line=~/dr(\d)_(\d+)/;
    $oo="HLA_DRB".$1."-".$2.".txt";
    open(TT, "$line") or die "$line";
    open(OUT, ">$oo") or die "$oo";
    print OUT "NumCols:\t9\n";
    while($ll=<TT>){
	@arr=split /\s+/, $ll;
	print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t0\t$arr[6]\t$arr[7]\t0\t$arr[9]\n";
    }
    print OUT "C\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}
