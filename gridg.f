      subroutine gridg ()

       implicit none
       include 'epm.inc'
       integer ig1,ig2,ig3,n1,n2,n3,ig

       ig=0
       do 10 ig3=0,ng3(is)-1
	  n3=ig3
	  if (ig3.gt.ng3(is)/2) n3=ig3-ng3(is)
	  do 20 ig2=0,ng2(is)-1
	     n2=ig2
	     if (ig2.gt.ng2(is)/2) n2=ig2-ng2(is)
	     do 30 ig1=0,ng1(is)-1
		n1=ig1
		if (ig1.gt.ng1(is)/2) n1=ig1-ng1(is)
		ig=ig+1
		pfft(n1,n2,n3)=ig
30           continue
20        continue
10     continue

       ngrid(is)=ig

       return
       end
