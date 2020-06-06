#! /usr/bin/perl

use strict;
use warnings;

sub read_file($$);
sub seqs_to_guide($$$);
sub check_structure($$$);
sub genetic_hash();
sub delete_gaps($$);
sub create_image($$$$);

# If there are no arguments, then the script fails: says how to use it and exit
if ($#ARGV == -1){
  print STDERR "ERROR: No rna file with sequences provided\n";
    print STDERR "Usage:\n ./mappin_rna.pl <input.txt> <output.txt> [VARNA.jar]\n";
  exit(1);
}
elsif ($#ARGV == 0){
  print STDERR "ERROR: No output file provided\n";
  print STDERR "Usage:\n ./mappin_rna.pl <input.txt> <output.txt> [VARNA.jar]\n";
  exit(1);
}

# If there's more one or more argument, uses the first one as the file with the RNA. The second one is a flag (-f) which tells if write in individual files or print the output
elsif ($#ARGV >= 0){
  my $file = $ARGV[0];
  my $output = $ARGV[1];
  my $varna = 0;
  my $varna_file_name;
  # Several ifs checks the number of arguments and if folders exists and conditions are meet
  # If file does not exists, exists throwing an error
  if (!(-e $file)){ print STDERR "File $file with RNA sequences can't be found\n"; exit(1);}
  
  # If there are 3 arguments, the last one is assumed to be the
  if ($#ARGV == 2){
      $varna_file_name = $ARGV[-1];
      $varna = 1;
      # The java file that creates the images
      # Checks if Varna file exists. If it does, puts a $varna flag on. Else, says that cannot be found and won't generate images
      if (!(-e $varna_file_name)){
        $varna = 0;
        print "WARNING: No VARNA.jar executable found, images won't be generated\n";
      }
    }
  # If theres no Varna file provided, it won't generate the images
  else {
    $varna = 0;
    }
  
  # Too many arguments is an unexpected behaviour, so throws an error and leaves.
  if ($#ARGV > 3){
    print STDERR "ERROR: Too many parameters indicated\n";
    print STDERR "Usage:\n ./mappin_rna.pl <input.txt> <output.txt> [VARNA.jar]\n";
    exit(1);
  }

  # Params are the name of the file with the sequences and the structure, which is the first element
  my $folder = "structure_images";
  # Images will be saved in $folder folder, and it will be created if it does not exists
  if (! (-e $folder) && $varna){mkdir($folder);};
  # $varna checks if images will (1) or won't (0) be created
  if ($varna){
    print "$varna_file_name executable found, images will be generated in \"$folder\" folder\n";
  };
    # Reading the sequences and saving them in %seqs
    my %seqs;
    print "Sequences and structures will be save on \"$output\" file\n";
    read_file($file, \%seqs);
    # Obtaining the main sequence and the reference structure in a hash
    my $refseq = $seqs{'SQref'};
    my $ref_linkage = $seqs{'SQstr'};
    # However, we will be iterating for the sequences so let's delete the structure and now, 
    # all values have the same structure (nucleotides)
    delete($seqs{'SQstr'});

    my %seqs_structures;
    my %guide;
    # Creating the structure guide that will help to get new structures
    seqs_to_guide($refseq, $ref_linkage, \%guide);

    open(OUT, ">$output");
    # For every sequence, let's calculate structure. clean it and save it
    foreach my $key (sort keys %seqs){
      print "Calculating structure of: $key\n";
      # Structure as an array and joining into a string
      my @seq_structure = check_structure($seqs{$key}, $ref_linkage, \%guide);
      my $structure = join("", @seq_structure);
      # Cleaning and deleting gaps
      delete_gaps(\$seqs{$key}, \$structure);
      # Writing into the file
      print OUT $key."\n".$seqs{$key}."\n".$structure."\n";
      if ($varna == 1){ # Checks if it has to create the image and do it if necessary
        print "Creating image of: $key\n";
        create_image($key, $seqs{$key}, $structure, $folder);
        };
      };
    # And the we create the images for them
};

# Read the file with the sequences
sub read_file($$){
  my $file = $_[0];
  my $seqs = $_[1];
  open(RNAs, $file);
  while (<RNAs>){
    my @line = split " ",  $_;
    $seqs->{$line[0]} = lc($line[1]);
  }
  close(RNAs);
}

# Turn the refseq structure to a dictionary that relates joined positions
sub seqs_to_guide($$$){
  my $refseq = $_[0];
  my $linkage = $_[1];
  my $guide = $_[2];
  my $depth = 0;
  my @depth_pos;

  # Foreach base, check if bonds and calculate depth, saving into the guide
  foreach (my $i = 0; $i < length($refseq); $i++){
    my $base_link = substr($linkage, $i, 1);
    if ($base_link eq "."){
      next;
      }
    elsif($base_link eq '(') {
      $depth_pos[$depth] = $i;
      $depth += 1
      }
    elsif ($base_link eq ')'){
      $depth -= 1;
      my $pos = "$depth_pos[$depth]";
      $$guide{$pos} = $i;
    };
    # Exit if finds an unpaired closing parenthesis
    if ($depth < 0){
      print STDERR "The reference structure is not correct in position $i\nClose parenthesis found without an equivelent open parenthesis";
      exit(-1);
    };
  };
  # Here exit if finds unpaired open parenthesis
  if ($depth > 0){
      print STDERR "The reference structure is not correct. At least $depth parenthesis are not closed";
      exit(-1);};
};

# Turn the reference structure to the other sequences
sub check_structure($$$){
  # Seq, Linkage array
  my $seq = $_[0];
  my @structure = split("", $_[1]);
  my $guide = $_[2];
  my @seq_structure;
  my %genetic_equivalence = &genetic_hash();
  for (my $i = 0; $i < length($seq); $i++){
    my $base = substr($seq, $i, 1);
    # Using the guide we get the paired base
    # If one base is lost, we should eliminate both its ( and the ) which pairs with
    if (($base eq '-') && ($structure[$i] eq '(')){
      push(@seq_structure, '.');
      # Buscar pareja y eliminarlo
      my $pos = $guide->{"$i"};
      $structure[$pos] = '.';
    }
    # Check if the pair still exists. If so, keep ')'. Else remove
    elsif ($structure[$i] eq '('){
      # Check here before pushing
      my $pair = $guide->{"$i"};
      my $base_pair = substr($seq, $pair, 1);
      if ($base_pair ne $genetic_equivalence{$base}){
        $structure["$i"] = '.';
        $structure["$pair"] = '.';
      }
      push(@seq_structure, $structure[$i]);

    }
    # If we find and ')', that means the equivalent '(' exists. In the previous elsif we took care of it
    # Technichally does the same as the else, but this way the idea is better understood. However, they could be fused to only one block
    elsif ($structure[$i] eq ')'){
      push(@seq_structure, $structure[$i]);
    }
    else{
        push(@seq_structure, $structure[$i]);
      }
  }
  return @seq_structure;
}

sub genetic_hash(){
  my %hash;
  $hash{'a'} = 'u';
  $hash{'u'} = 'a';
  $hash{'c'} = 'g';
  $hash{'g'} = 'c';
  return %hash;
}

sub delete_gaps($$){
  my $seq = $_[0];
  my $structure = $_[1];
  my @new_seq;
  my @new_struct;
  # IF it's a dash, delete the base and associated structure.
  for(my $i = 0; $i < length($$seq); $i++){
    if (substr($$seq, $i, 1) ne '-'){
      push(@new_seq, substr($$seq, $i, 1));
      push(@new_struct, substr($$structure, $i, 1));}
  }
  $$seq = join("", @new_seq);
  $$structure = join("", @new_struct);
}

sub create_image($$$$){
  my $ID = $_[0];
  my $seq = $_[1];
  my $struct = $_[2];
  my $folder = $_[3];
  my $file = $folder."/".$ID.".png";
  system("java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd  -sequenceDBN $seq -structureDBN \"$struct\" -o $file 2>> Varna.log >> Varna.log");
}
