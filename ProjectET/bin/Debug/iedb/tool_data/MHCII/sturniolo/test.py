from bisect import bisect, bisect_left, bisect_right
a=open("perc_HLA_DRB1-0101.txt", 'r').readlines()

ll=[]
for  i in range(len(a)):
	ll.append(float(a[i].strip()))

tt=-12.5

print ll[0], ll[9999]
print bisect(ll, tt)
print bisect_left(ll, tt)
print ll.index(-12.5)
