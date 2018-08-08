SUBROUTINE handle_command(outbuffer, buffer, buffer_length) BIND ( C )

  USE, INTRINSIC :: ISO_C_BINDING
  USE RUN
  IMPLICIT NONE

  CHARACTER(KIND=C_CHAR, LEN=1), INTENT(IN)  :: buffer(1024)
  INTEGER(C_INT), INTENT(IN) :: buffer_length
  CHARACTER(KIND=C_CHAR, LEN=1), INTENT(OUT) :: outbuffer(1024)

!  INTERACTIVE_MODE=.TRUE.
!  INTERUPT=.TRUE.

  print *,"You said '",buffer(1:buffer_length),"' at TIME = ", TIME

!  WRITE(outbuffer,"(1(A))") "handled command."

!  WRITE(outchar,"(3(A),F10.5,1(A))")"You said '",buffer(1:buffer_length),"' at TIME = ", TIME, CHAR(0)

END SUBROUTINE handle_command
