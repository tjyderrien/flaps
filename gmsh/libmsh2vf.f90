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
MODULE LIBMSH2VF
	IMPLICIT NONE
	CONTAINS
!
	SUBROUTINE print_help()
		print*, "msh2vf options"
		print*, "--------------"
		print*, "-i : precise the input file (Type of = string)"
		print*, "-o : precise the output file (Type of input = string) - OPTIONAL"
		print*, " "
		print*, "Examples:    msh2vf -i my_mesh.msh"
		print*, "             msh2vf -i my_mesh.msh -o my_mesh.vf"
		stop
	END SUBROUTINE print_help
!
!
	SUBROUTINE extract_parameters(namefile_msh, namefile_vf)
		CHARACTER(LEN=*), INTENT(INOUT) :: namefile_msh, namefile_vf
!
		INTEGER          :: nbarg, ierror, i
		CHARACTER(LEN=2) :: flag
		LOGICAL          :: input_given, output_given
!
		nbarg = iargc()
		input_given = .FALSE.
		output_given = .FALSE.
		i = 1
		CALL getarg(i, flag)
		i = i+1
		DO WHILE (i <= nbarg)
			IF (flag(2:2) == "i") THEN
				! Specify input
				CALL getarg(i, namefile_msh)
				input_given = .TRUE.
				i = i+1
			ELSE IF (flag(2:2) == "o") THEN
				! Specify output
				CALL getarg(i, namefile_vf)
				output_given = .TRUE.
				i = i+1
			END IF
			IF (i < nbarg) THEN
				CALL getarg(i, flag)
				i = i+1
			END IF
		END DO
		IF (.NOT. (input_given)) THEN
			CALL print_help()
		ELSE
			IF (.NOT. (output_given)) THEN
				namefile_vf = namefile_msh(1:(LEN(TRIM(namefile_msh))-4)) // ".vf"
				print*, "No output given. Output file will be ", TRIM(namefile_vf)
			END IF
		END IF
	END SUBROUTINE extract_parameters
!
!
	RECURSIVE SUBROUTINE sort_dicho_int(tab)
		INTEGER, DIMENSION(:,:), INTENT(INOUT) :: tab
!
		INTEGER                               :: i, j, k, n, ierror
		INTEGER, DIMENSION(1:1,1:size(tab,2)) :: temp
		INTEGER, DIMENSION(:,:), POINTER      :: tab_temp
!   
		n = size(tab,1)
		IF (n == 2) THEN
			IF (tab(2,1) < tab(1,1)) THEN
				temp(1,:) = tab(1,:)
				tab(1,:) = tab(2,:)
				tab(2,:) = temp(1,:)
			ELSE IF ((tab(2,1) == tab(1,1)) .AND. (tab(2,2) < tab(1,2))) THEN
				temp(1,:) = tab(1,:)
				tab(1,:) = tab(2,:)
				tab(2,:) = temp(1,:)
			END IF
		ELSE IF (n > 2) THEN
			CALL sort_dicho_int(tab(1:int(n/2),:))
			CALL sort_dicho_int(tab(int(n/2)+1:n,:))
       		allocate(tab_temp(size(tab,1),size(tab,2)), STAT=ierror)
			tab_temp = tab
			i = 1
			j = int(n/2) + 1
			k = 1
			DO WHILE ((i <= int(n/2)) .OR. (j <= n))
				IF (i > int(n/2)) THEN
					tab(k,:) = tab_temp(j,:)
					j=j+1
					k=k+1
				ELSE IF (j > n) THEN
					tab(k,:) = tab_temp(i,:)
					i=i+1
					k=k+1
				ELSE IF (tab_temp(i,1) < tab_temp(j,1)) THEN
					tab(k,:) = tab_temp(i,:)
					i = i+1
					k = k+1
				ELSE IF ((tab_temp(i,1) == tab_temp(j,1)) .AND. (tab_temp(i,2) < tab_temp(j,2))) THEN
					tab(k,:) = tab_temp(i,:)
					i = i+1
					k = k+1
				ELSE
					tab(k,:) = tab_temp(j,:)
					j = j+1
					k = k+1
				END IF
			END DO
			deallocate(tab_temp)
		END IF
	END SUBROUTINE sort_dicho_int
!
!
	SUBROUTINE correct_orientation(vertices, cells)
		DOUBLE PRECISION, DIMENSION(:,:), POINTER, INTENT(IN) :: vertices
		INTEGER, DIMENSION(:,:), POINTER, INTENT(INOUT)       :: cells
!
		INTEGER          :: i, itemp
		DOUBLE PRECISION :: xA, xB, xC, yA, yB, yC, det1, det2
!
		DO i=1,size(cells,1)
			IF (size(cells,2) == 4) THEN
				! We treat triangles
				xA = vertices(cells(i,1),1)
				xB = vertices(cells(i,2),1)
				xC = vertices(cells(i,3),1)
				yA = vertices(cells(i,1),2)
				yB = vertices(cells(i,2),2)
				yC = vertices(cells(i,3),2)
				det1 = xA*yB-xB*yA + xB*yC-xC*yB + xC*yA-xA*yC
				IF (det1 < 0d0) THEN
					! The triangle is counter-trigonometric oriented -> Need to correct it
					itemp = cells(i,1)
					cells(i,1) = cells(i,2)
					cells(i,2) = itemp
					!print*, "Needed to change the orientation of triangle ", i
				END IF
			ELSE IF (size(cells,2) == 5) THEN
				xA = vertices(cells(i,1),1)
				xB = vertices(cells(i,2),1)
				xC = vertices(cells(i,3),1)
				yA = vertices(cells(i,1),2)
				yB = vertices(cells(i,2),2)
				yC = vertices(cells(i,3),2)
				det1 = xA*yB-xB*yA + xB*yC-xC*yB + xC*yA-xA*yC
				xA = vertices(cells(i,2),1)
				xB = vertices(cells(i,3),1)
				xC = vertices(cells(i,4),1)
				yA = vertices(cells(i,2),2)
				yB = vertices(cells(i,3),2)
				yC = vertices(cells(i,4),2)
				det2 = xA*yB-xB*yA + xB*yC-xC*yB + xC*yA-xA*yC
				IF ((det1 < 0d0) .AND. (det2 < 0d0)) THEN
					! The quadrangle is counter-trigonometric oriented -> Need to correct it
					itemp = cells(i,1)
					cells(i,1) = cells(i,4)
					cells(i,4) = itemp
					itemp = cells(i,2)
					cells(i,2) = cells(i,3)
					cells(i,3) = itemp
					!print*, "Needed to change the orientation of quadrangle ", i
				ELSE IF (det1*det2 < 0d0) THEN
					print*, "Quadrangle ", i, " may be not well built"
					print*, "The program will stop now..."
					stop
				END IF
			ELSE
				print*, "Only triangles or quadrangles can be treated by SUBROUTINE correct_orientation."
				print*, "The program will stop now..."
				stop
			END IF
		END DO
	END SUBROUTINE correct_orientation
!
!
	SUBROUTINE read_msh_file(namefile, vertices, points, segments, triangles, quadrangles, dim_physical_entities, &
					& id_physical_entities, name_physical_entities, idvertices)
		DOUBLE PRECISION, DIMENSION(:,:), POINTER, INTENT(INOUT) :: vertices
		INTEGER, DIMENSION(:), POINTER, INTENT(OUT)		   :: idvertices
		INTEGER, DIMENSION(:,:), POINTER, INTENT(INOUT)          :: points, segments, triangles, quadrangles
		INTEGER, DIMENSION(:), POINTER, INTENT(INOUT)            :: dim_physical_entities, id_physical_entities
		CHARACTER(LEN=*), DIMENSION(:), POINTER, INTENT(INOUT)   :: name_physical_entities
		CHARACTER(LEN=*), INTENT(IN)                             :: namefile
!
		CHARACTER(LEN=300)               :: buf
		INTEGER                          :: nb_physical_entities
		INTEGER                          :: nb_elements
		INTEGER                          :: nb_vertices
		INTEGER                          :: nb_points, ip
		INTEGER                          :: nb_segments, is
		INTEGER                          :: nb_triangles, it
		INTEGER                          :: nb_quadrangles, iq
		INTEGER                          :: ierror, i, itemp1, itemp2, itemp3
		INTEGER, DIMENSION(:), POINTER   :: element_type
		INTEGER, DIMENSION(:,:), POINTER :: element_tag
		DOUBLE PRECISION                 :: dtemp1
!
		OPEN(UNIT=1, FILE=TRIM(namefile), STATUS="unknown")
		! Read $MeshFormat (3 lines)
		READ(UNIT=1,FMT=*)
		READ(UNIT=1,FMT=*)
		READ(UNIT=1,FMT=*)
		! Read $PhysicalNames if present, or $Nodes
		READ(UNIT=1,FMT=*) buf
		IF (TRIM(buf) == "$PhysicalNames") THEN
			! Read the number of declared physical entities (> 0)
			READ(UNIT=1,FMT=*) nb_physical_entities
			IF (nb_physical_entities < 1) THEN
				print*, "The section $PhysicalNames should contains at least one entry"
				stop
			END IF
			allocate(dim_physical_entities(1:nb_physical_entities), STAT=ierror)
			allocate(id_physical_entities(1:nb_physical_entities), STAT=ierror)
			allocate(name_physical_entities(1:nb_physical_entities), STAT=ierror)
			! Read the list of declared physical entities
			DO i=1,nb_physical_entities
				READ(UNIT=1,FMT=*) dim_physical_entities(i), id_physical_entities(i), name_physical_entities(i)
			END DO
			! Read $EndPhysicalNames
			READ(UNIT=1,FMT=*)
			! Read $Nodes
			READ(UNIT=1,FMT=*)
		END IF
		! Read the number of vertices
		READ(UNIT=1,FMT=*) nb_vertices
		! Read the list of vertices
		allocate(vertices(1:nb_vertices,1:3), STAT=ierror)
		allocate(idvertices(1:nb_vertices), STAT=ierror)
		DO i=1,nb_vertices
			READ(UNIT=1, FMT=*) idvertices(i), vertices(i,1), vertices(i,2), dtemp1
		END DO
		! Read $EndNodes
		READ(UNIT=1,FMT=*)
		! Read $Elements
		READ(UNIT=1,FMT=*)
		! Read the number of elements
		READ(UNIT=1,FMT=*) nb_elements
		! Read the type and the tags for each element
		nb_points = 0
		nb_segments = 0
		nb_triangles = 0
		nb_quadrangles = 0
		allocate(element_type(1:nb_elements), STAT=ierror)
		allocate(element_tag(1:nb_elements,4), STAT=ierror)
		DO i=1,nb_elements
			READ(UNIT=1,FMT=*) itemp1, element_type(i), element_tag(i,1)
			IF (element_tag(i,1) < 2) THEN
				print*, "Each element requires at least a physical tag and a geometry tag. Please correct the input .msh or .geo files."
				stop
			ELSE IF (element_tag(i,3) == 3) THEN
				print*, "Mesh partitioning is not handled. Please correct the input .msh or .geo files."
				stop
			END IF
			SELECT CASE (element_type(i))
				CASE (1, 8, 26, 27, 28)
					! The element is a segment
					nb_segments = nb_segments+1
				CASE (2, 9, 20, 21, 22, 23, 24, 25)
					! The element is a triangle
					nb_triangles = nb_triangles+1
				CASE (3, 10, 16)
					! The element is a quadrangle
					nb_quadrangles = nb_quadrangles+1
				CASE (15)
					! The element is a point
					nb_points = nb_points+1
				CASE DEFAULT
					print*, "The input file contains some 3D elements. It cannot be handled."
					stop
			END SELECT
		END DO
		! Go back onto the beginning of the elements list
		CLOSE(UNIT=1)						! Close the file
		OPEN(UNIT=1, FILE=TRIM(namefile), STATUS="unknown")	! Open the file
		READ(UNIT=1,FMT=*)					! Read $MeshFormat
		READ(UNIT=1,FMT=*)					!
		READ(UNIT=1,FMT=*)					! Read $EndMeshFormat
		READ(UNIT=1,FMT=*) buf					! Read $PhysicalNames or $Nodes
		IF (TRIM(buf) == "$PhysicalNames") THEN
			READ(UNIT=1,FMT=*)				! Read the number of physical entities
			DO i=1,nb_physical_entities			! Read the list of physical entities
				READ(UNIT=1,FMT=*)			! |
			END DO						! |_
			READ(UNIT=1,FMT=*)				! Read $EndPhysicalNames
			READ(UNIT=1,FMT=*)				! Read $Nodes
		END IF
		READ(UNIT=1,FMT=*)					! Read the number of nodes
		DO i=1,nb_vertices					! Read the list of nodes
			READ(UNIT=1, FMT=*)				! |
		END DO							! |_
		READ(UNIT=1,FMT=*)					! Read $EndNodes
		READ(UNIT=1,FMT=*)					! Read $Elements
		READ(UNIT=1,FMT=*)					! Read the number of elements
		IF (nb_points > 0) THEN
			allocate(points(1:nb_points, 1:2), STAT=ierror)
		END IF
		IF (nb_segments > 0) THEN
			allocate(segments(1:nb_segments, 1:3), STAT=ierror)
		END IF
		IF (nb_triangles > 0) THEN
			allocate(triangles(1:nb_triangles, 1:4), STAT=ierror)
		END IF
		IF (nb_quadrangles > 0) THEN
			allocate(quadrangles(1:nb_quadrangles, 1:5), STAT=ierror)
		END IF
		ip = 1
		is = 1
		it = 1
		iq = 1
		DO i=1,nb_elements
			SELECT CASE (element_type(i))
				CASE (15)
					! The element is a point
					READ(UNIT=1,FMT=*) itemp1, itemp2, itemp3, element_tag(i,2:(element_tag(i,1)+1)), points(ip,1)
					points(ip,2) = element_tag(i,2)
					ip = ip+1
				CASE (1, 8, 26, 27, 28)
					! The element is a segment
					READ(UNIT=1,FMT=*) itemp1, itemp2, itemp3, element_tag(i,2:(element_tag(i,1)+1)), segments(is,1:2)
					segments(is,3) = element_tag(i,2)
					is = is+1
				CASE (2, 9, 20, 21, 22, 23, 24, 25)
					! The element is a triangle
					READ(UNIT=1,FMT=*) itemp1, itemp2, itemp3, element_tag(i,2:(element_tag(i,1)+1)), triangles(it,1:3)
					triangles(it,4) = element_tag(i,2)
					it = it+1
				CASE (3, 10, 16)
					! The element is a quadrangle
					READ(UNIT=1,FMT=*) itemp1, itemp2, itemp3, element_tag(i,2:(element_tag(i,1)+1)), quadrangles(iq,1:4)
					quadrangles(iq,5) = element_tag(i,2)
					iq = iq+1
			END SELECT
		END DO
		CLOSE(UNIT=1)
		deallocate(element_type)
		deallocate(element_tag)
	END SUBROUTINE read_msh_file
!
!
	SUBROUTINE compute_edges(triangles, quadrangles, edges)
		INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)          :: triangles, quadrangles
		INTEGER, DIMENSION(:,:), POINTER, INTENT(INOUT)       :: edges
!
		INTEGER                                   :: i, ierror, n1, n2, k, nbtriangles, nbquadrangles, nbedges
		INTEGER, DIMENSION(:,:), POINTER          :: temp1, temp2, tempedges
		DOUBLE PRECISION, DIMENSION(:,:), POINTER :: temp3, temp4
		IF (associated(triangles)) THEN
			nbtriangles = size(triangles,1)
		ELSE
			nbtriangles = 0
		END IF
		IF (associated(quadrangles)) THEN
			nbquadrangles = size(quadrangles,1)
		ELSE
			nbquadrangles = 0
		END IF
		nbedges = 3*nbtriangles+4*nbquadrangles
		allocate(tempedges(1:nbedges,1:12), STAT=ierror)
		tempedges(:,:) = 0
		n1 = 0
		n2 = 0
		! We consider two temporary edges lists: temp1 will contains the edges with vertex indices (i1,i2) satisfying i1 < i2 and
		! temp2 will contains the edges with vertex indices (i1,i2) satisfying i1 > i2
		! Each edge can appear only one time in temp1
		! Each edge can appear only one time in temp2
		! An interior edge appears one time in temp1 and one time in temp2
		! A bound edge appears one time in the reunion of temp1 and temp2
		! In temp2, we permute i1 and i2 before storing the edge in order to have temp2(*,1) < temp2(*,2)
		! ==> temp1 and temp2 are not completely fullfilled after having read all the cells
		allocate(temp1(1:nbedges,1:5), STAT=ierror)
		allocate(temp2(1:nbedges,1:5), STAT=ierror)
		temp1(:,:) = 0
		temp2(:,:) = 0
		DO i=1,nbtriangles
			IF (triangles(i,1) < triangles(i,2)) THEN
				temp1(n1+1,1) = triangles(i,1)	! First extremity of the edge
				temp1(n1+1,2) = triangles(i,2)	! Second extremity of the edge
				temp1(n1+1,3) = i		! Index of involved triangular cell
				temp1(n1+1,4) = triangles(i,3)	! Keep in mind the third vertex of the triangle
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = triangles(i,2)	! First extremity of the edge
				temp2(n2+1,2) = triangles(i,1)	! Second extremity of the edge
				temp2(n2+1,3) = i		! Index of involved triangular cell
				temp2(n2+1,4) = triangles(i,3)	! Keep in mind the third vertex of the triangle
				n2 = n2+1
			END IF
			IF (triangles(i,2) < triangles(i,3)) THEN
				temp1(n1+1,1) = triangles(i,2)	! First extremity of the edge
				temp1(n1+1,2) = triangles(i,3)	! Second extremity of the edge
				temp1(n1+1,3) = i		! Index of involved triangular cell
				temp1(n1+1,4) = triangles(i,1)	! Keep in mind the third vertex of the triangle
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = triangles(i,3)	! First extremity of the edge
				temp2(n2+1,2) = triangles(i,2)	! Second extremity of the edge
				temp2(n2+1,3) = i		! Index of involved triangular cell
				temp2(n2+1,4) = triangles(i,1)	! Keep in mind the third vertex of the triangle
				n2 = n2+1
			END IF
			IF (triangles(i,3) < triangles(i,1)) THEN
				temp1(n1+1,1) = triangles(i,3)	! First extremity of the edge
				temp1(n1+1,2) = triangles(i,1)	! Second extremity of the edge
				temp1(n1+1,3) = i		! Index of involved triangular cell
				temp1(n1+1,4) = triangles(i,2)	! Keep in mind the third vertex of the triangle
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = triangles(i,1)	! First extremity of the edge
				temp2(n2+1,2) = triangles(i,3)	! Second extremity of the edge
				temp2(n2+1,3) = i		! Index of involved triangular cell
				temp2(n2+1,4) = triangles(i,2)	! Keep in mind the third vertex of the triangle
				n2 = n2+1
			END IF
		END DO
		! In the case of triangles, we duplicate the 4th column of temp1 and temp2
		temp1(1:n1,5) = temp1(1:n1,4)
		temp2(1:n2,5) = temp2(1:n2,4)
		DO i=1,nbquadrangles
   	  		IF (quadrangles(i,1) < quadrangles(i,2)) THEN
				temp1(n1+1,1) = quadrangles(i,1)	! First extremity of the edge
				temp1(n1+1,2) = quadrangles(i,2)	! Second extremity of the edge
				temp1(n1+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp1(n1+1,4) = quadrangles(i,4)	! Keep in mind the vertex which preceeds temp1(:,1)
				temp1(n1+1,5) = quadrangles(i,3)	! Keep in mind the vertex which follows temp1(:,2)
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = quadrangles(i,2)	! First extremity of the edge
				temp2(n2+1,2) = quadrangles(i,1)	! Second extremity of the edge
				temp2(n2+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp2(n2+1,4) = quadrangles(i,3)	! Keep in mind the vertex which follows temp2(:,1)
				temp2(n2+1,5) = quadrangles(i,4)	! Keep in mind the vertex which preceeds temp2(:,2)
				n2 = n2+1
			END IF
			IF (quadrangles(i,2) < quadrangles(i,3)) THEN
				temp1(n1+1,1) = quadrangles(i,2)	! First extremity of the edge
				temp1(n1+1,2) = quadrangles(i,3)	! Second extremity of the edge
				temp1(n1+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp1(n1+1,4) = quadrangles(i,1)	! Keep in mind the vertex which preceeds temp1(:,1)
				temp1(n1+1,5) = quadrangles(i,4)	! Keep in mind the vertex which follows temp1(:,2)
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = quadrangles(i,3)	! First extremity of the edge
				temp2(n2+1,2) = quadrangles(i,2)	! Second extremity of the edge
				temp2(n2+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp2(n2+1,4) = quadrangles(i,4)	! Keep in mind the vertex which follows temp2(:,1)
				temp2(n2+1,5) = quadrangles(i,1)	! Keep in mind the vertex which preceeds temp2(:,2)
				n2 = n2+1
			END IF
			IF (quadrangles(i,3) < quadrangles(i,4)) THEN
				temp1(n1+1,1) = quadrangles(i,3)	! First extremity of the edge
				temp1(n1+1,2) = quadrangles(i,4)	! Second extremity of the edge
				temp1(n1+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp1(n1+1,4) = quadrangles(i,2)	! Keep in mind the vertex which preceeds temp1(:,1)
				temp1(n1+1,5) = quadrangles(i,1)	! Keep in mind the vertex which follows temp1(:,2)
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = quadrangles(i,4)	! First extremity of the edge
				temp2(n2+1,2) = quadrangles(i,3)	! Second extremity of the edge
				temp2(n2+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp2(n2+1,4) = quadrangles(i,1)	! Keep in mind the vertex which follows temp2(:,1)
				temp2(n2+1,5) = quadrangles(i,2)	! Keep in mind the vertex which preceeds temp2(:,2)
				n2 = n2+1
			END IF
			IF (quadrangles(i,4) < quadrangles(i,1)) THEN
				temp1(n1+1,1) = quadrangles(i,4)	! First extremity of the edge
				temp1(n1+1,2) = quadrangles(i,1)	! Second extremity of the edge
				temp1(n1+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp1(n1+1,4) = quadrangles(i,3)	! Keep in mind the vertex which preceeds temp1(:,1)
				temp1(n1+1,5) = quadrangles(i,2)	! Keep in mind the vertex which follows temp1(:,2)
				n1 = n1+1
			ELSE 
				temp2(n2+1,1) = quadrangles(i,1)	! First extremity of the edge
				temp2(n2+1,2) = quadrangles(i,4)	! Second extremity of the edge
				temp2(n2+1,3) = nbtriangles+i		! Index of involved quadrangular cell
				temp2(n2+1,4) = quadrangles(i,2)	! Keep in mind the vertex which follows temp2(:,1)
				temp2(n2+1,5) = quadrangles(i,3)	! Keep in mind the vertex which preceeds temp2(:,2)
				n2 = n2+1
			END IF
		END DO
		! We sort temp1 and temp2
		CALL sort_dicho_int(temp1)
		CALL sort_dicho_int(temp2)
		! We compute the number of false edges in temp1 and temp2
		n1 = 1
		DO WHILE (temp1(n1,1) == 0)
			n1 = n1+1
		END DO
		n2 = 1
		DO WHILE (temp2(n2,2) == 0)
			n2 = n2+1
		END DO
		!==> The n1 first lines of temp1 and the n2 first lines of temp2 are made of 0. We skip these lines.
		! We fill the 8 first columns of edges
		i = 1
		DO WHILE ((n1 <= nbedges) .OR. (n2 <= nbedges))
			IF (n1 > nbedges) THEN
				! temp1 has been completely scanned
				! The data entirely comes from temp2
				tempedges(i,1:3) = temp2(n2,1:3)
				tempedges(i,5:6) = temp2(n2,4:5)
				n2 = n2+1
				i = i+1
			ELSE IF (n2 > nbedges) THEN
				! temp2 has been completely scanned
				! The data entirely comes from temp1 
				tempedges(i,1:3) = temp1(n1,1:3)
				tempedges(i,5:6) = temp1(n1,4:5)
				n1 = n1+1
				i = i+1
			ELSE IF ((temp1(n1,1) == temp2(n2,1)) .AND. (temp1(n1,2) == temp2(n2,2))) THEN
				! the edges from temp1 and temp2 coincide --> we treat an interior edge
				! The extremities come from temp1
				tempedges(i,1:2) = temp1(n1,1:2)
				! The index of the cell on the left side (K) of the edge comes from temp1
				tempedges(i,3) = temp1(n1,3)
				! The index of the cell on the right side (L) of the edge comes from temp2
				tempedges(i,4) = temp2(n2,3)
				! The indices of preceeding and following vertices in K come from temp1
				tempedges(i,5:6) = temp1(n1,4:5)
				! The indices of preceeding and following vertices in L come from temp1
				tempedges(i,7:8) = temp2(n2,4:5)
				i = i+1
				n1 = n1+1
				n2 = n2+1
			ELSE IF ((temp1(n1,1) == temp2(n2,1)) .AND. (temp1(n1,2) < temp2(n2,2))) THEN
				! The line from temp1 preceeds the line from temp2 in the sorting convention
				! --> The data entirely comes from temp1
				tempedges(i,1:3) = temp1(n1,1:3)
				tempedges(i,5:6) = temp1(n1,4:5)
				i = i+1
				n1 = n1+1
			ELSE IF (temp1(n1,1) < temp2(n2,1)) THEN
				! The line from temp1 preceeds the line from temp2 in the sorting convention
				! --> The data entirely comes from temp1
				tempedges(i,1:3) = temp1(n1,1:3)
				tempedges(i,5:6) = temp1(n1,4:5)
				i = i+1
				n1 = n1+1
			ELSE
				! The line from temp2 preceeds the line from temp1 in the sorting convention
				! --> The data entirely comes from temp2
				tempedges(i,1:3) = temp2(n2,1:3)
				tempedges(i,5:6) = temp2(n2,4:5)
				i = i+1
				n2 = n2+1
			END IF
		END DO
		! Extract the real list of edges
		nbedges = 1
		DO WHILE (.NOT.((tempedges(nbedges,1) == 0) .AND. (tempedges(nbedges,2) == 0)))
			nbedges = nbedges+1
		END DO
		nbedges = nbedges-1
		allocate(edges(1:nbedges,1:13), STAT=ierror)
		edges(1:nbedges,1:12) = tempedges(1:nbedges,1:12)
		! edges(i,13) will contain the identifier of the physical zone in which the edge is included (1D or 2D subdomain).
		! At present time, it is fixed to 0.
		edges(1:nbedges,13) = 0
		deallocate(temp1)
		deallocate(temp2)
		deallocate(tempedges)
	END SUBROUTINE compute_edges
!
!
	SUBROUTINE identify_physical_zone_for_edges(edges, segments, triangles, quadrangles)
		INTEGER, DIMENSION(:,:), POINTER, INTENT(INOUT) :: edges
		INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)    :: segments, triangles, quadrangles
!
		INTEGER                          :: i, ierror, itemp, is
		INTEGER                          :: nb_segments, nb_triangles, nb_quadrangles, nb_edges
		INTEGER                          :: K, L, K_phys_id, L_phys_id
		INTEGER, DIMENSION(:,:), POINTER :: copy_segments
!
		allocate(copy_segments(1:size(segments,1),1:size(segments,2)), STAT=ierror)
		copy_segments(:,:) = segments(:,:)
		CALL sort_dicho_int(copy_segments)
		DO i=1,size(copy_segments,1)
			IF (copy_segments(i,1) > copy_segments(i,2)) THEN
				itemp = copy_segments(i,1)
				copy_segments(i,1) = copy_segments(i,2)
				copy_segments(i,2) = itemp
			END IF
		END DO
		CALL sort_dicho_int(copy_segments)
		IF (associated(triangles)) THEN
			nb_triangles = size(triangles,1)
		ELSE
			nb_triangles = 0
		END IF
		IF (associated(quadrangles)) THEN
			nb_quadrangles = size(quadrangles,1)
		ELSE
			nb_quadrangles = 0
		END IF
		is = 1
		DO i=1,size(edges,1)
			K = edges(i,3)
			L = edges(i,4)
			IF (L > 0) THEN
				! We treat an interior edge
				IF (K > nb_triangles) THEN
					K_phys_id = quadrangles(K-nb_triangles,5)
				ELSE
					K_phys_id = triangles(K,4)
				END IF
				IF (L > nb_triangles) THEN
					L_phys_id = quadrangles(L-nb_triangles,5)
				ELSE
					L_phys_id = triangles(L,4)
				END IF			
				IF (K_phys_id == L_phys_id) THEN
					! K and L are included in the same physical 2D subdomain
					edges(i,13) = K_phys_id
				ELSE
					! K and L are in different physical subdomains which are separated by a physical 1D subdomain
					! in which edge(i,:) is included
					! edges(i,1:2) should appear in the segments list
					DO WHILE (.NOT.((copy_segments(is,1) == edges(i,1)).AND.(copy_segments(is,2) == edges(i,2))))
						is = is+1
					END DO
					edges(i,13) = copy_segments(is,3)
				END IF
			ELSE
				! We treat a boundary edge
				! edges(i,1:2) should appear in the segments list
				DO WHILE (.NOT.((copy_segments(is,1) == edges(i,1)).AND.(copy_segments(is,2) == edges(i,2))))
					is = is+1
				END DO
				edges(i,13) = copy_segments(is,3)
				edges(i,4) = -copy_segments(is,3)
			END IF
		END DO
		deallocate(copy_segments)
	END SUBROUTINE identify_physical_zone_for_edges
!
!
	SUBROUTINE extract_boundedges(edges, boundedges)
		INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)    :: edges
		INTEGER, DIMENSION(:,:), POINTER, INTENT(INOUT) :: boundedges
!
		INTEGER :: nbboundedges, i, j, ierror
!
		nbboundedges = 0
		DO i=1,size(edges,1)
			IF (edges(i,4) <= 0) THEN
				nbboundedges = nbboundedges+1
			END IF
		END DO
		allocate(boundedges(1:nbboundedges,1:4), STAT=ierror)
		j = 1
		DO i=1,size(edges,1)
			IF (edges(i,4) <= 0) THEN
				boundedges(j,1:2) = edges(i,1:2)
				boundedges(j,3) = edges(i,4)
				boundedges(j,4) = i
				j = j+1
			END IF
		END DO
	END SUBROUTINE extract_boundedges
!
!
	SUBROUTINE save_vf(vertices, triangles, quadrangles, boundedges, edges, namefile_vf)
		DOUBLE PRECISION, DIMENSION(:,:), POINTER, INTENT(IN) :: vertices
		INTEGER, DIMENSION(:,:), POINTER, INTENT(IN)          :: triangles, quadrangles, boundedges, edges
		CHARACTER(LEN=*), INTENT(IN)                          :: namefile_vf
!
		INTEGER                           :: i, ierror, ios
		INTEGER                           :: nb_vertices, nb_triangles, nb_quadrangles, nb_edges, nb_boundedges
		CHARACTER(LEN=LEN(namefile_vf)-3) :: namemesh
!
		namemesh(1:LEN(namemesh)) = namefile_vf(1:LEN(namemesh))
		OPEN(UNIT=1, FILE=namefile_vf, FORM="formatted", ACCESS="sequential", STATUS="replace", ACTION="write", &
					& POSITION="rewind", IOSTAT=ios)
		nb_vertices = size(vertices,1)
		IF (associated(triangles)) THEN
			nb_triangles = size(triangles,1)
		ELSE
			nb_triangles = 0
		END IF
		IF (associated(quadrangles)) THEN
			nb_quadrangles = size(quadrangles,1)
		ELSE
			nb_quadrangles = 0
		END IF
		nb_edges = size(edges,1)
		nb_boundedges = size(boundedges,1)
! 		WRITE(*,FMT=*) "Maillage cree avec msh2vf"
! 		WRITE(*,FMT=*) " "
! 		WRITE(*,FMT=*) " "
! 		WRITE(*,FMT=*) "Version"
! 		WRITE(*,FMT=*) "4"
! 		WRITE(*,FMT=*) "Nom du maillage"
! 		WRITE(*,FMT=*) TRIM(namemesh)
		WRITE(*,FMT=*) "MESH INFORMATIONS"
! 		WRITE(*,FMT=*) "Nombre de triangles"
! 		WRITE(*,FMT=*) nb_triangles
		WRITE(*,FMT=*) "Quadrangles number"
		WRITE(*,FMT=*) nb_quadrangles
		WRITE(*,FMT=*) "Edges"
		WRITE(*,FMT=*) nb_edges
		WRITE(*,FMT=*) "Boundedges"
		WRITE(*,FMT=*) nb_boundedges
		WRITE(*,FMT=*) "Nombre de sommets"
		WRITE(*,FMT=*) nb_vertices
		WRITE(*,FMT=*) "Choix des centres"
		WRITE(*,FMT=*) .FALSE.
		WRITE(*,FMT=*) "Sommets"
		DO i=1,nb_vertices
			WRITE(*,FMT=*) vertices(i,1:2)
		END DO
! 		WRITE(UNIT=1,FMT=*) "Arretes"
! 		DO i=1,nb_edges
! 			WRITE(UNIT=1,FMT=*) edges(i,1:8), 0, 0, 0, 0
! 		END DO
! 		WRITE(UNIT=1,FMT=*) "Triangles"
! 		DO i=1,nb_triangles
! 			WRITE(UNIT=1,FMT=*) triangles(i,1:3)
! 		END DO
! 		WRITE(UNIT=1,FMT=*) "Quadrangles"
! 		DO i=1,nb_quadrangles
! 			WRITE(UNIT=1,FMT=*) quadrangles(i,1:4)
! 		END DO
! 		WRITE(UNIT=1,FMT=*) "Arretes du bord"
! 		DO i=1,nb_boundedges
! 			WRITE(UNIT=1,FMT=*) boundedges(i,1:4)
! 		END DO
! 		CLOSE(UNIT=1)
	END SUBROUTINE save_vf
!
!
	SUBROUTINE convert_to_vf(namefile_msh, namefile_vf)
		CHARACTER(LEN=*), INTENT(IN)               :: namefile_msh, namefile_vf
!
		DOUBLE PRECISION, DIMENSION(:,:), POINTER :: vertices
		INTEGER, DIMENSION(:,:), POINTER          :: points, segments, triangles, quadrangles, boundedges, edges
		INTEGER, DIMENSION(:), POINTER            :: dim_physical_entities, id_physical_entities, idvertices
		CHARACTER(LEN=200), DIMENSION(:), POINTER :: name_physical_entities
		INTEGER                                   :: k
		CHARACTER(LEN=1)                          :: choice
!	
		! We read the .msh file
		CALL read_msh_file(namefile_msh, vertices, points, segments, triangles, quadrangles, dim_physical_entities, &
					& id_physical_entities, name_physical_entities, idvertices)
		! We verify if the triangles are sorted in trigonometric sense and we bring correction if necessary
		IF (associated(triangles)) THEN
			CALL correct_orientation(vertices, triangles)
		END IF
		! We verify if the quadrangles are sorted in trigonometric sense and we bring correction if necessary
		IF (associated(quadrangles)) THEN
			CALL correct_orientation(vertices, quadrangles)
		END IF
		! We build the edges
		CALL compute_edges(triangles, quadrangles, edges)
		! For each edge, we identify the physical zone in which it is included
		CALL identify_physical_zone_for_edges(edges, segments, triangles, quadrangles)
		! We extract the bound edges
		CALL extract_boundedges(edges, boundedges)
		! We save the .vf file
! 		CALL save_vf(vertices, triangles, quadrangles, boundedges, edges, namefile_vf)
! 		print*, "File ", TRIM(namefile_vf), " saved."
	END SUBROUTINE convert_to_vf
END MODULE LIBMSH2VF
