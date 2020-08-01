
from ip import IP	
import time

start_time = time.time()
angles = [0,15,30,45,60] #vinkler
angles = [60]

for i in range(len(angles)):

	c=IP(Ncol=260, Nrow=344, theta=angles[i])
	c.invasion()
	print("--- %s seconds ---" % (time.time() - start_time))