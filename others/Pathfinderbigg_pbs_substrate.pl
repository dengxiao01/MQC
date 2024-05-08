#FILE: Pathfinderbigg_pbs.pl
#AUTH: Xiaoping Liao (liao_xp@tib.cas.cn)
#DATE: Dec. 1st, 2022
#VERS: 1.0

# This pipeline is prepared to automatically do pathfinder for many products.
use warnings;
use strict;
use File::Temp;
use File::Path;
use File::Basename;
use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

my $command="perl Pathfinderbigg_pbs_iJN1463.pl @ARGV";

my %options = (
    substrate                  => '/public/home/liaoxp/Projects/Pathfinder/substrates.txt',
    product                    => '/public/home/liaoxp/Projects/Pathfinder/products.txt',
    scripts                     => '/public/home/liaoxp/Projects/Pathfinder/',
    verbose                     => undef,
    help                        => undef
);

GetOptions(
    'sub=s'         => \$options{substrate},
    'scripts=s'    => \$options{scripts},
    'pro=s'          => \$options{product},
    'o=s'          => \$options{output},
    'v'            => \$options{verbose},
    'h|help'       => \$options{help}
);

if ( defined( $options{help} ) ) {
    print_usage();
    exit(0);
}

# check the options
check_option(
    option_name     => 'pro',
    option_value    => $options{product},
    must_be_defined => 1
);
check_option(
    option_name     => 'sub',
    option_value    => $options{substrate},
    must_be_defined => 1
);
check_option(
    option_name     => 'o',
    option_value    => $options{output},
    must_be_defined => 1
);
check_option(
    option_name     => 'scripts',
    option_value    => $options{scripts},
    must_be_defined => 0
);




# check the directories
mkdir $options{output} unless -d $options{output};

mkpath("$options{output}/Results") unless -d "$options{output}/Results";
mkpath("$options{output}/Results/Logs") unless -d "$options{output}/Results/Logs";
mkpath("$options{output}/Results/output_pbs") unless -d "$options{output}/Results/output_pbs";

message( $options{verbose},$options{log_file},"The full Command used is: $command\n" );


# $sample_number is used to count the number of samples
my $sample_number = 0;

open my $SUB,'<'.$options{substrate};
my @sub=();
while(<$SUB>){
    chomp;
    push @sub, $_;
}

open my $PRO,'<'.$options{product};
my @pro=();
while(<$PRO>){
    chomp;
    push @pro, $_;
}

my $job_count=0;
my $str='';
foreach  my $p (@pro){
    foreach  my $s (@sub){
        #print $s."\n";
           #print $m."\n";
           #$job_count++;
           $str=$s.",".$p;
           $str =~ s/\r//g;
           #print $str."\n";
           open my $PBS, '>'. "$options{output}/Results/output_pbs/$str.pbs";
           print $PBS '#PBS -N '.$str."\n";
           print $PBS '#PBS -o '."$options{output}".'/Results/Logs/'.$str.".log\n";
           print $PBS "#PBS -l walltime=10:00:00\n";
           print $PBS "#PBS -l nodes=1:ppn=1\n";
           print $PBS "#PBS -q high\n";
           print $PBS "#PBS -j oe\n";
           print $PBS "source activate\n";
           print $PBS "conda activate cplex\n";
           print $PBS "cd /public/home/wei_f/bigg_pathway\n";
           #print $PBS 'python pathfind.py models/'."\n";
           my $ss="$s $p";
           $ss =~ s/\r//g;
           $p =~ s/\r//g;
           print $ss."\t$options{output}\n";
           print $PBS "python pathfindbigg_substrate.py $ss $options{output}  > $options{output}".'/Results/'."$p.txt\n";
           print $PBS "conda deactivate\n";
           $job_count++;
           `qsub $options{output}/Results/output_pbs/$str.pbs`; 
    }
}
message( $options{verbose}, undef, "There are total $job_count Jobs in this queue and all the typical analysis have been submitted to the cluster. You can find the results in the output folder $options{output}/Results. Of course, it is not finished yet. For each sample, there is a output folder output_pbs. You can find all the .pbs that are submitted, also you can find the corresponding log file for each submmitted job.\n" );







sub check_option {
    my %args = (@_);

    if ( $args{must_be_defined} ) {
        if ( !defined( $args{option_value} ) ) {
            print "Option '$args{option_name}' must be defined.\n";
            print_usage();
            exit(1);
        }
    }

    if ( !defined( $args{option_value} ) ) {
        return;
    }

    if ( $args{type_int} ) {
        if ( $args{option_value} =~ m/[^\d\+\-]/ ) {
            die("Option '$args{option_name}' must be an integer.");
        }
    }

    if ( $args{type_real} ) {
        if ( $args{option_value} =~ m/[^\d\+\-\.]/ ) {
            die("Option '$args{option_name}' must be a real number.");
        }
    }

    if ( defined( $args{max_value} ) ) {
        if ( $args{option_value} > $args{max_value} ) {
            die("Option '$args{option_name}' must be less than or equal to $args{max_value}."
            );
        }
    }

    if ( defined( $args{min_value} ) ) {
        if ( $args{option_value} < $args{min_value} ) {
            die("Option '$args{option_name}' must be greater than or equal to $args{min_value}."
            );
        }
    }

    if ( defined( $args{allowed_values} ) ) {
        my $is_allowed = 0;
        foreach my $allowed_value ( @{ $args{allowed_values} } ) {
            if ( lc( $args{option_value} ) eq lc($allowed_value) ) {
                $is_allowed = 1;
                last;
            }
        }
        if ( !$is_allowed ) {
            die("Option '$args{option_name}' must match one of the allowed values: "
                    . join( ",", @{ $args{allowed_values} } )
                    . "." );
        }
    }
}



sub message {
    my $verbose = shift;
    my $log     = shift;
    my $message = shift;
    $message =~ s/\s(\s{1,})/$1/g;
    if ($verbose) {
        print $message;
    }
    if ($log) {
        print $log $message;
    }
}


sub print_usage {
    print <<BLOCK;
USAGE:
perl Pathfinderbigg_pbs.pl [-arguments]

 -sub [File]        : Input file containing substrates (Required).
 -pro [File]        : output file containing products (Required).
 -o [Dir]        : output folder (Required).
 -v              : Provide progress messages (Optional).
 -help           : Provide list of arguments and exit (Optional).
 
BLOCK
}
