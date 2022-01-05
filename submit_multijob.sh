bsub -P acc_als-omics -q long -n 4 -W 72:00 -R span[hosts=1] -R rusage[mem=3750] -L /bin/bash -oo multijob.stdout -eo multijob.stderr  "sh script_multijob.sh"
