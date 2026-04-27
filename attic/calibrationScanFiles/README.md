# Create Pulseq and TOPPE scan files for B0 mapping


## Dependencies

Requires the 
[Pulseq](http://pulseq.github.io/)
,
[TOPPE](https://toppemri.github.io/)
, and
[PulseGEq](https://github.com/toppeMRI/PulseGEq/)
packages.


Get the code:
```
$ cd ~/github/
$ git clone git@github.com:pulseq/pulseq.git
$ git clone git@github.com:toppeMRI/toppe.git
% git clone git@github.com:toppeMRI/PulseGEq.git
```

Set Matlab paths:
```
>> addpath ~/github/pulseq/matlab/            % +mr package
>> addpath ~/github/toppe/                    % +toppe package
>> addpath ~/github/PulseGEq/                 % +pulsegeq package
```


## Get the code and create scan files

```
$ cd ~/github/
$ git clone git@github.com:HarmonizedMRI/B0shimming.git
```

```
>> makescanfiles;
```




