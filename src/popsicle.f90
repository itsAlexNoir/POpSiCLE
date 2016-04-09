
 !     d8888b.  .d88b.  d8888b. .d8888. d888888b  .o88b. db      d88888b
 !     88  `8D .8P  Y8. 88  `8D 88'  YP   `88'   d8P  Y8 88      88'
 !     88oodD' 88    88 88oodD' `8bo.      88    8P      88      88ooooo
 !     88~~~   88    88 88~~~     `Y8b.    88    8b      88      88~~~~~
 !     88      `8b  d8' 88      db   8D   .88.   Y8b  d8 88booo. 88.
 !     88       `Y88P'  88      `8888Y' Y888888P  `Y88P' Y88888P Y88888P
 !

MODULE POPSICLE

  ! Here the list of modules of the library
  USE constants
  USE tools
  USE bessel
  USE gaussleg
  USE scatt_interp
  USE cubic_interp
  USE scattcoords
  USE cubcoords
  USE coords
  USE waveduplicate
  USE scattboundcyl
  USE scattboundcart
  USE cubboundcyl
  USE cubboundcart
  USE boundary
  USE fourier
  USE sht
  USE patchwork
  USE io_pop
  USE io_subset
  USE io_surface
  USE surface
  USE flux
  USE observables
  USE io_sampt

  IMPLICIT NONE


END MODULE POPSICLE
