      subroutine skip(ifile,ilines)
      integer ifile,ilines,i
      do i=1,ilines
         read(ifile,*)
      end do
      end
