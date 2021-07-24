#!/usr/bin/perl

# Regex to match any floating point number
$fpr="[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?";

# List of velocities to scan over
@V=(0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.15,1.35);
@V=(1,2,3,4,6,8,10,12,14,16,18);

# List of viscosities
@nu=(2.5,4,6.5,10,13,16,20,25,32,40,50,65,100,160);

# Create a sequence of configuration files
foreach $k (0..$#nu) {

    # Switch domain size
    $fn=$k<7?"lo":"ld";

    # Open master configuration file and output file
    open A,"template.cfg" or die "Can't open master configuration file\n";
    open B,">${fn}_sig0_vis$nu[$k].cfg" or die "Can't open output file\n";

    # Search and replace
    while(<A>) {

        # Match keyword, spaces, and a floating point number. Replace with the
        # given velocity.
        s/^(nul_cSt +)$fpr/$1$nu[$k]/;
        #s/^(V +)$fpr/${1}1.05/;

        # For the larger viscosities, increase the simulation size and duration
        if($k>=7) {
            s/^(L_nd +)$fpr/${1}45/;
            s/^(t_end_nd +)$fpr/${1}100/;
        }

        print B;
    }
    close B;
    close A;
}
