import os
import time

for i in range( 1, 16 ): 
	a = str(1000 + (i - 1) * 400 + 1)
	b = str(1000 + (i)*400 )
	c = "sbatch --output=/dev/null --error=/dev/null --array=" + a + "-" + b + " ccbar.slurm.sh"
	print c
	os.system( c )
	time.sleep( 60 * 60 )
