
ECOCLIMAP data is downloaded from ECOCLIMAP web-page:
http://www.cnrm.meteo.fr/gmme/PROJETS/ECOCLIMAP/page_ecoclimap.htm

The original data consists of these files/directories:
------------------------------------------------------
convert_lai.sh
DATA
ECOCLIMAP_LICENCE_AGREEMENT.doc
README.txt
script_compile
sources

Here the sources directory includes the original source code for 
the data base according to "october 04, corrected september 2006".

Additional files/directories added at Rossby Centre:
----------------------------------------------------
README_Rossby_Centre
RUN
SCRIPTS

The RUN directory includes modified code compared to the sources directory:
---------------------------------------------------------------------------
- Redefinition of NTILE=3 in vegtype_to_patch.f90
- Added output of F_lake in ecoclimap_parameters.f90
- Reading and processing of soil carbon in ecoclimap_parameters.f90

The results of processed ECOCLIMAP data are stored in RUN/GLOBAL_NSCALE/
where NSCALE refers to the resolution (see ecoclimap.f90)

The variables stored belongs to four categories:
------------------------------------------------
- No dependence on time and tile
- Dependence on tile
	t01	open land low vegetation
	t02	coniferous forest
	t03	broadleaf forest
- Dependence on time for month=1-12
- Dependence on time and tile 

Here the variables per category are listed:

No dependence on time and tile
------------------------------
albedo_soil.NSCALE
clay.NSCALE
F_bare.NSCALE
f_broadleaf_tree.NSCALE
f_C3_crops.NSCALE
f_C4_crops.NSCALE
f_conif_tree.NSCALE
f_flat_bare_soil.NSCALE
f_grassland.NSCALE
f_irr_crops.NSCALE
F_lowv.NSCALE
f_park_marshes.NSCALE
f_perm_snow.NSCALE
frac_tile01.NSCALE
frac_tile02.NSCALE
frac_tile03.NSCALE
f_rocks.NSCALE
F_tree.NSCALE
f_trop_grass.NSCALE
f_tropical_tree.NSCALE
F_urb.NSCALE
F_wat.NSCALE
sand.NSCALE
F_lake.NSCALE
soil_carb.NSCALE

Dependence on tile
------------------
albedo_veg.tile.NSCALE
d_root.tile.NSCALE
d_soil.tile.NSCALE
rsmin.tile.NSCALE

Dependence on time
------------------
glob_alb.time.NSCALE

Dependence on time and tile
---------------------------
alb.time.tile.NSCALE
emis.time.tile.NSCALE
greeness.time.tile.NSCALE
lai.time.tile.NSCALE
veg.time.tile.NSCALE
z0.time.tile.NSCALE
z0h.time.tile.NSCALE

----------------------------------------------
----------------------------------------------

ECOCLIMAP as used by RCA:
-------------------------

In GEMINI.f the array eco(KLON,KLAT,meco) is used for ECOCLIMAP info
as used in the physics package (SL2TIM and below). Only a subset
of the variables above are used in RCA. Currently meco=43 and
refers to the following lists of variables. The variable names are
those as used in PHYS.f

The first part (nn=1-9) concerns those variables that vary with time of the year.
Thus, these variables should be continiously updated during a simulation using linear
interpolation between monthly values.
														Missing data in
c     variable  description                                     nn	ECOCLIMAP file				ECOCLIMAP file (see "Missing data" below!)
c     --------  -----------------------------------------       --	--------------				------------------------------------------
c     lai_t1	  LAI open land					 1	time_of_year( lai.time.t01.NSCALE )	999
c     lai_t2	  LAI conif forest				 2	time_of_year( lai.time.t02.NSCALE )	999
c     lai_t3	  LAI broad-leaf forest				 3	time_of_year( lai.time.t03.NSCALE )	999
c     z0_t1	  roughness length open land			 4	time_of_year( z0.time.t01.NSCALE )	999  See "Roughness" below!
c     z0_t2	  roughness length conif forest			 5	time_of_year( z0.time.t02.NSCALE )	999  See "Roughness" below!
c     z0_t3	  roughness length broad-leaf forest		 6	time_of_year( z0.time.t03.NSCALE )	999  See "Roughness" below!
c     emis_t1	  emissivity open land				 7	time_of_year( emis.time.t01.NSCALE )	999
c     alb_t1	  albedo open land				 8	time_of_year( alb.time.t01.NSCALE )	999
c     veg_t1	  vegetation cover open land			 9	time_of_year( veg.time.t01.NSCALE )	999

The second part (nn=10-meco) concerns those variables that are constant with time.
														Missing data in
c     variable  description                                     nn	ECOCLIMAP file				ECOCLIMAP file (see "Missing data" below!)
c     --------  -----------------------------------------       --	--------------				------------------------------------------
c     frland	  fraction land					10	1 - F_wat.NSCALE			-
c     alb_soil	  albedo of bare soil				11	albedo_soil.NSCALE			999
c     clay	  percentage of clay				12	clay.NSCALE				where F_wat.NSCALE =1
c     sand	  percentage of sand				13	sand.NSCALE				where F_wat.NSCALE =1
c     frac_t1	  fraction open land				14	frac_tile01.NSCALE			999
c     frac_t2	  fraction conif forest				15	frac_tile02.NSCALE			999
c     frac_t3	  fraction broad-leaf forest			16	frac_tile03.NSCALE			999
c     alb_t2	  albedo conif forest				17	alb.07.t02.NSCALE			999
c     alb_t3	  albedo broad-leaf forest			18	alb.07.t03.NSCALE			999
c     emis_t2	  emissivity conif forest			19	emis.07.t02.NSCALE			999
c     emis_t3	  emissivity broad-leaf forest			20	emis.07.t03.NSCALE			999
c     veg_t2	  vegetation cover conif forest			21	veg.07.t02.NSCALE			999
c     veg_t3	  vegetation cover broad-leaf forest		22	veg.07.t03.NSCALE			999
c     droot_t1	  root depth open land				23	d_root.t01.NSCALE			999
c     droot_t2	  root depth conif forest			24	d_root.t02.NSCALE			999
c     droot_t3	  root depth broad-leaf forest			25	d_root.t03.NSCALE			999
c     dsoil_t1	  soil depth open land				26	d_soil.t01.NSCALE			999
c     dsoil_t2	  soil depth conif forest			27	d_soil.t02.NSCALE			999
c     dsoil_t3	  soil depth broad-leaf forest			28	d_soil.t03.NSCALE			999
c     rsmin_t1	  minimum surface resistance open land		29	rsmin.t01.NSCALE			999
c     rsmin_t2	  minimum surface resistance conif forest	30	rsmin.t02.NSCALE			999
c     rsmin_t3 	  minimum surface resistance			31	rsmin.t03.NSCALE			999
c     alb_veg_t1  albedo vegetation open land			32	albedo_veg.t01.NSCALE			999
c     alb_veg_t2  albedo vegetation conif forest		33	albedo_veg.t02.NSCALE			999
c     alb_veg_t3  albedo vegetation				34	albedo_veg.t03.NSCALE			999
c     texture	  texture according to texture triangle		35  Processed info see "Texture" below		-
c     minlai_t1	  annual min of LAI open land			36	min( lai.time.t01.NSCALE )		999
c     minlai_t2	  annual min of LAI conif forest		37	min( lai.time.t02.NSCALE )		999
c     minlai_t3	  annual min of LAI broad-leaf forest		38	min( lai.time.t03.NSCALE )		999
c     maxlai_t1	  annual max of LAI open land			39	max( lai.time.t01.NSCALE )		999
c     maxlai_t2	  annual max of LAI conif forest		40	max( lai.time.t02.NSCALE )		999
c     maxlai_t3	  annual max of LAI broad-leaf forest		41	max( lai.time.t03.NSCALE )		999
c     frac_lake	  fraction of lake				42	F_lake.NSCALE				-
c     soil_carb   soil carbon                                   43 	soil_carb.NSCALE			where F_wat.NSCALE =1


Missing data
------------
All grid boxes with missing data should be excluded in aggregation.

Roughness
---------
The roughness length should be aggregated using this method:

Texture
-------
The texture is calculated as a function sand and clay.
The texture should be calculated on aggregated sand and clay information (sand.NSCALE and clay.NSCALE).
Two functions are used to classify a soil in the texture triangle based on percentage of sand and clay:
gimle: /home/sm_aulle/ECOCLIMAP/w.f  and /home/sm_aulle/ECOCLIMAP/i.f


