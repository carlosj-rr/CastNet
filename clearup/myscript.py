import sys

args=sys.argv
if len(args) > 1:
	run_prefix=sys.argv[1]+"-"
	print(f"Run prefix provided: {run_prefix}")
else:
	print("A run prefix was not provided, using default \'RUN-\'")
	run_prefix="RUN-"

print('run prefix:',run_prefix)
