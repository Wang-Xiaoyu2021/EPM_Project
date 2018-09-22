c     This subroutine reads in a line from the specified file.
c     It ignores all lines beginning with the comment characters # and !
c     and moves onto the next line.
      
c     ierr returns the status of the attempt to read a line
c     0 => successful reading of a non-comment line
c     1 => end of file reached
c     2 => blank line
      
      subroutine getline(ifile,linestring,ierr)
      implicit none
      
      integer ifile,ierr
c      character*72 linestring
      character*120 linestring
      logical iread
      
      ierr=0
      iread=.true.
      do while (iread)
         read (ifile,'(a)',end=100) linestring
            if (linestring(1:1)/='#' .and. linestring(1:1)/='!') then
               iread=.false.
               if (linestring(1:1)==' ') ierr=2
            end if
      end do
      
      return
      
 100  ierr=1
      
      end
