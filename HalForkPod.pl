#!/usr/bin/perl
#
# cactus batch run 
#


use Parallel::ForkManager;
my $max = 4;
my $pm = Parallel::ForkManager->new($max);


my %ids = ('TpodH1' => 't_pod_1',
	'TpodE240140H1' => 't_pod_24_0140_h1',
	'TpodE240140H2' => 't_pod_24_0140_h2',
	'TpodE240154H1' => 't_pod_24_0154_h1',
	'TpodE240154H2' => 't_pod_24_0154_h2');

FILES:
foreach $file (@ARGV){
	$pm->start and next FILES; ## fork
	$file =~ m/cactusPodura_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)/ or die "failed here $file\n";
	$id1 = $ids{$1};
	$id2 = $ids{$2};
	$out = "out_syn_$1"."_$2".".psl";
	system "~/source/hal/bin/halSynteny --queryGenome $id2 --targetGenome $id1 $file $out\n";
	
	$pm->finish;
}

$pm->wait_all_children;



