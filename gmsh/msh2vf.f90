!----------------------------------------------------------------------------!
! Copyright (c) 2014 Alexandre MOUTON alexandre.mouton@math.univ-lille1.fr   !
!                                                                            !
! This file is part of MESH2VF.                                              !
!                                                                            !
! MESH2VF is free software: you can redistribute it and/or modify it under   !
! the terms of the GNU General Public License as published by the Free       !
! Software Foundation, either version 3 of the License, or at your option)   !
! any later version.                                                         !
!                                                                            !
! MESH2VF is distributed in the hope that it will be useful, but WITHOUT ANY !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or          !
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for   !
! more details.                                                              !
!                                                                            !
! You should have received a copy of the GNU General Public License along    !
! with MESH2VF. If not, see <http://www.gnu.org/licenses/>.                  !
!----------------------------------------------------------------------------!
PROGRAM MSH2VF
	USE LIBMSH2VF
!
	CHARACTER(LEN=300) :: namefile_vf, namefile_msh
!
	CALL extract_parameters(namefile_msh, namefile_vf)
	CALL convert_to_vf(namefile_msh, namefile_vf)
END PROGRAM MSH2VF
