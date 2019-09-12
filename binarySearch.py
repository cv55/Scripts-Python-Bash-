
##Function that searches for elements in huge lists combining binary search and recursion. 
##The list must be sorted to apply binary search.


def binarySearch(lis,l,h,tofind):
	lis=sorted(lis)
	m=int((l+h)/2)
	if l==h:
		if lis[l]!=tofind:
			return "Value {0} not found in list {1}".format(tofind,lis)
		else:
			return "The value {0} is found with index {1} in the sorted list {2}".format(tofind, l, lis)
	if lis[m]>tofind:
		h=m-1
	elif lis[m]<tofind:
		l=m+1
	else:
		return "The value {0} is found with index {1} in the sorted list {2}".format(tofind, m, lis)

	return binarySearch(lis,l,h,tofind)




a=[22,33,11,6,9,100,99]
print(binarySearch(a,len(a)-len(a),len(a)-1,98))
