import sys
import re
import operator


fn = sys.argv[1]
#p = re.compile("TEST")
test_count = {}
with open(fn) as f:
    for line in f:
        match_group = re.match(r'TEST (\d+) (\d+)', line, re.M|re.I)
        test_num = match_group.group(1)
        count = match_group.group(2)
        #print"{}: {}".format(test_num, count)
        test_count[test_num] = int(count)

#for key, value in test_count.iteritems():
    #print value
sorted_names = sorted(test_count, key=lambda x:test_count[x])
for k in sorted_names:
    #print("{} : {}".format(k, test_count[k]))
    #print("%d "%format(k), end='')
    sys.stdout.write('%s '%k)
    
    sys.stdout.flush()

