#!/usr/bin/perl
#
# Linear and high-order shim calibration scan on GE MRI scanners. 
# 
# Prescribe the TOPPE sequence in ../psd/Cartesian/, do Auto-prescan (APS) and Download, 
# and then run this script (from anywhere on the scanner host):
# 	$ perl shimcal.pl
#
# Jon-Fredrik Nielsen 2019, University of Michigan

# system(". /usr/local/bin/portableCmds_sh");

# acquisition loop parameters
@Shims = ('x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy');
@AmpLinear = (-10, 10);       # loop over these (integer) hardware amplitudes for each linear shim
@AmpHO = (-500, 500);         # loop over these (integer) hardware amplitudes for each higher-order shim

$atpfile = 'tmp.atp';

# initialize shim amplitude array
@ShimAmps = @Shims;
for ($i = 0; $i < @ShimAmps; $i++) {
	@ShimAmps[$i] = 0;
}

# scan 
for ($i = 0; $i < @Shims; $i++) {

	if ($i > 0) {
	#	 last;    # exit the for loop
	}

	print	"Setting $Shims[$i] shim and scanning.\n";

	# get appropriate amplitude array for this shim term
	if ( @Shims[$i] eq 'x' || @Shims[$i] eq 'y' || @Shims[$i] eq 'z' ) {
		@Amp=@AmpLinear;
	} else {
		@Amp=@AmpHO;
	}

	# For each shim amp: set shims and run one scan (acquire one 3D B0 map)
	for($j = 0; $j < @Amp; $j++) {

		@ShimAmps[$i] = @Amp[$j];

		# linear shim setting command string
		$slin = '';
		for($ii = 0; $ii < 3; $ii++) {
			$slin .= "\"@ShimAmps[$ii]\"";
			if ($ii < 2) {
				$slin .= " ";
			}
		}

		# high-order shim setting command string
		$sho = '';
		$shoZero = '';
		for($ii = 3; $ii < @Shims; $ii++) {
			$sho .= "@Shims[$ii] @ShimAmps[$ii] ";
			$shoZero .= "@Shims[$ii] 0 ";
		}

		print	"\tLinear shims:$slin\n";
		print "\tHO shims: $sho\n";

		# run the scan
		#system("readallcurrent");
		system("/usr/g/bin/setNavShimCurrent $sho");
		sleep(5);
		# atpexec("SYSTEM \"/usr/g/bin/setNavShimCurrent $sho\"");
		# setCurrentLinearShim @ShimAmps[0] @ShimAmps[1] @ShimAmps[2]
		atpexec("PSC_SET_SHIM $slin");
		sleep(15);     # NB! Seems to be important to have some delay here for the shim values to 'take'
		atpexec("ACCEPT_RX");
		atpexec("DOWNLOAD");
		### atpexec("APS_EVENT");
		atpexec("SCAN_EVENT");
		atpexec("RECON_STOPPED");
		sleep(2);                            # allow some time for scanner to start writing Pfile
		system("/usr/g/bin/serviceClassA OFF"); 

		# Wait for Pfile to be written, then rename it.
		print "\tWaiting for pfile...";

		$datdir = "/usr/g/mrraw";

		$pfileSize = 0;
		$lastPfileSize = -1;
		while ($pfileSize - $lastPfileSize > 0) {
			$lastPfileSize = $pfileSize;
			sleep(2); 
			@Tmp = `ls -tr $datdir/P[0-9][0-9][0-9][0-9][0-9].7 | tail -1`;
			#chomp(@Tmp);         # remove newline character at end of string
			@Tmp = `du -b @Tmp[0]`;
			($pfileSize, $pfile) = split(/\s+/, @Tmp[0]);  # delimiter \s+ = any number of white spaces
		}

		$newName = "$datdir/P,@Shims[$i],@Amp[$j].7";
		system("mv $pfile $newName "); 
		print " done\n";
	}

	# Set amplitude for this shim term back to zero
	@ShimAmps[$i] = 0;   
}

system("/usr/g/bin/setNavShimCurrent $shoZero");
system("/usr/g/bin/serviceClassA OFF"); 

# TODO: Move Pfiles to server 


# execute one ATP command
sub atpexec {
	$cmd = $_[0];
	$atpfile = "tmp.atp";
	open(fh, '>', "$atpfile");
	print fh "$cmd;\n";
	close(fh);
	system("cat $atpfile");
	system("/usr/g/bin/serviceClassA ON"); 
	sleep(1);
	system("atp $atpfile");
	sleep(1);
}
