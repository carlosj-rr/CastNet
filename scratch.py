import retikulos as rt
import cProfile, pstats, io
np=rt.np
from pstats import SortKey
founder=rt.founder_miner()
pr = cProfile.Profile()
pr.enable()
new_pop=rt.grow_pop(founder,100)
pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
textfile=open("profile_growpop_new2.txt","w")
textfile.write(s.getvalue())
textfile.close()