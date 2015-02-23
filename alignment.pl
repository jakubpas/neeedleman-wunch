#!/usr/bin/perl
# (c) 2001 by Jakub Pas.
# Perl Implementation of Neeedleman-Wunch Affine Gap Algorithm

if ($#ARGV!=1) {die "Usage alignment.pl seq1 seq2";}

my $gapopen=1; # koszt wstwienia gapu
my $gapext=1; # koszt rozszerzenia gapu
my $inf=100;
my @x; # sekwencja 1
my @y; # sekwencja 2
my @xl; #aligned x
my @yl; #aligned y
my @score; # Score matrix
my @direction; # Traceback matrix
my @matrix; # substitution matrix
my @matrow;

#affine gap
my @scoregapx;
my @scoregapy;

&readmatrix("/cygdrive/c/Program Files/BioEdit/tables/BLOSUM62", @matrix,@matrow);

# &drawtable (@matrow, @matrow, @matrix);

unshift (@matx,"");
unshift (@maty,"");

my $seq1=&fastatoarray($ARGV[0],@x);
my $seq2=&fastatoarray($ARGV[1],@y);

unshift (@x,"");
unshift (@y,"");

# Inicjalizacja tablic;

$score[0][0] = 0;
$scoregapx[0][0]= $scoregapy[0][0] = 0;
for ($i=1;$i<=$#x;$i++) {$scoregapy[$i][0] = &w($i);}
for ($j=1;$j<=$#y;$j++) {$scoregapy[0][$j] = &w($j);}
for ($j=1;$j<=$#y;$j++) {$scoregapx[0][$j] = $inf;}
for ($i=1;$i<=$#x;$i++) {$scoregapy[$i][0] = $inf;}

# Policz score dla matrycy

for ($i=1;$i<=$#x;$i++) {
    for ($j=1;$j<=$#y;$j++) {
        $scoregapx[$i][$j] = &min( $score[$i-1][$j] + &w(1), $scoregapx[$i-1][$j] + $gapopen );
        $scoregapy[$i][$j] = &min( $score[$i][$j-1] + &w(1), $scoregapy[$i][$j-1] + $gapopen );
        $score[$i][$j] = &min( $score[$i-1][$j-1] + &getscore($x[$i],$y[$j],@matrix), &min($scoregapx[$i][$j], $scoregapy[$i][$j]));
    }
    
    
}
#    &drawtable (@x, @y, @score);
#    &drawtable (@x, @y, @scoregapx);
#    &drawtable (@x, @y, @scoregapy);


# Drukuj tablice;

#&drawtable (@x, @y, @direction);
#&drawtable (@x, @y, @score);

# Traceback

#&traceback (@score,@direction,$i-1,$j-1,@xl,@yl);
&traceback2 (@score,@scoregapx,@scoregapy,$i,$j,@xl,@yl);

@xl=reverse(@xl);
@yl=reverse(@yl);

print ">".$seq1."n";
print @xl;
print "n";
print ">".$seq2."n";
print @yl;
print "n";

exit;

#subs;

sub traceback {
    
    my ($points,$kierunek,$i,$j,$xl,$yl) = @_;
    
    if ($$kierunek[$i][$j]==1) {push (@xl, $x[$i]); push (@yl,$y[$j]);$i--;$j--; &traceback (@score,@direction,$i,$j,@xl,@yl);}
    elsif ($$kierunek[$i][$j]==2) {push (@xl, $x[$i]);push (@yl, "-"); $i--; &traceback (@score,@direction,$i,$j,@xl,@yl);}
    elsif ($$kierunek[$i][$j]==3) {push (@yl, $y[$j]);push (@xl, "-"); $j--; &traceback (@score,@direction,$i,$j,@xl,@yl);}
    
}

sub filetostring {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /.gz|.Z/) {
        if (! open (FILE, "gzip -dc $file |")) { print "no file $file";}
    }
    elsif (! open (FILE, $file)) { print  "no file $filen";}
    local $buf = <FILE>;
    close (FILE);
    return $buf;
}

sub filetoarray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /.gz|.Z/) {
        if (! open (FILE, "gzip -dc $file |")) { print "no file $file";}
    }
    elsif (! open (FILE, $file)) { print  "no file $filen";}
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf=split(/n/,$buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

sub maximum {
    my $max;
    
    
    return $max;
}

sub drawtable {
    
    my ($x, $y, $table) = @_;
    
    my $i,$j;
    
    print "   ";
    for ($i=0;$i<=@$x-1;$i++) { print $$x[$i]."|";}
    print "n";
    
    for ($j=0;$j<=@$y-1;$j++) {
        for ($i=0;$i<=@$x-1;$i++)
        {
            if ($i==0 and $j==0) {print " ";}
            if ($i==0) {print $$y[$j]."|";}
            print $$table[$i][$j]."|";
        }
        print "n";
    }
    print "n";
}

sub cleartable {
    
    my ($lx, $ly, $table) = @_;
    my $i, $j;
    
    for ($j=0;$j<=$ly;$j++) {
        for ($i=0;$i<=$lx;$i++) {
            $$table[$i][$j]=0;
        }
    }
    
}

sub fastatoarray {
    
    my($filename,$seq)=@_;
    my $tmpseq;
    
    $tmpseq=&filetostring($filename);
    $tmpseq=~s/>(.*)n(.*)/$2/;
    my $name=$1;
    $tmpseq=~s/[a-z]/[A-Z]/g;
    $tmpseq=~s/[^A-Z]//g;
    @seq=split('',$tmpseq);
    my $x=0;
    foreach $_ (@seq) {$$seq[$x]=$_; $x++}
    return $name;
}

sub readmatrix {
    
    ($filename, $matrix, $matrow) = @_;
    my $y=0;
    @mtxfile=&filetoarray("$filename");
    foreach $_ (@mtxfile) {
        @mtxrow = split (" ",$_);
        if ($#mtxrow == 23) {
            @$matrow=@mtxrow;
        }
        if ($#mtxrow == 24) {
            shift (@mtxrow);
            for ($x=0;$x<=23;$x++) {
                $$matrix[$x][$y]=$mtxrow[$x];
            }
            $y++;
        }
    }
}

sub getscore {
    
    ($tmpx,$tmpy,$matrix) = @_;
    my $x,$y;
    my $i=0;
    foreach $_ (@matrow) {
        if ($tmpx eq $_) {$x=$i}
    $i++}
    my $i=0;
    foreach $_ (@matrow) {
        if ($tmpy eq $_) {$y=$i}
    $i++}
    return $$matrix[$x][$y];
    
}

sub max {
    
    my ($x, $y) = @_;
    
    if ($x>=$y) { return $x;}
    else {return $y;}
    
}

sub min {
    
    my ($x, $y) = @_;
    
    if ($x<=$y) { return $x;}
    else {return $y;}
    
}

sub w {
    
    my $i=shift;
    
    return $gapopen + ($gapext * $i)
    
}

sub traceback2 {
    
    my ($points,$left,$up,$i,$j,$xl,$yl) = @_;
    
    
    while ( $i != 1 and $j != 1 )
    {
        if ( $$points[$i][$j] == $$left[$i][$j])
        {
            $i--;
            push (@xl, $x[$i]); push (@yl, "-");
        }
        elsif ( $$points[$i][$j] == $$points[$i-1][$j-1] + &getscore($i,$j,@matrix) )
        {
            $i--; $j--;
            push (@xl, $x[$i]); push (@yl, $x[$i]);
        }
        elsif ( $$points[$i][$j] == $$up[$i][$j])
        {
            $j--;
            push (@yl, $y[$i]); push (@xl, "-");
        }
    }
    while ( $i != 1 ) {
        $i--;
        push (@xl, $x[$i]); push (@yl, "-");
    }
    while ( $j != 1 ) {
        $j--;
        push (@yl, $y[$i]); push (@xl, "-")
    }
}