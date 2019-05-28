from bisect import insort
lists=open("HLA_list.txt", 'r').readlines()
for i in range(len(lists)):
	arr=[]
	ff="swissprot_"+lists[i].strip()
	scan = open(ff, 'r').readlines()
	for j in range(len(scan)):
		temp=scan[j].strip()
		insort(arr, temp)
	print arr[0], arr[1]
