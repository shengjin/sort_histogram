
ledebut := $(shell echo '')
labranche := $(shell basename ${PWD})
lemilieu:= $(shell echo '_')
lafin := $(shell svn info | grep "Last Changed Rev" | cut -f 2 -d: | sed -e 's/^[ \t]*//')
#lafin :=$(shell echo '21')
codename := ${ledebut}${labranche}${lemilieu}${lafin}

hostname := $(shell hostname -s)
#standart compiler
compiler := gfortran
warning := true
#warning := false

DF=
DF90=
DL=
# special options for some hosts
ifeq "$(hostname)" "isis"
   compiler := ifort8.0
   FC=/opt/intel_fc_80/bin/ifort
   F90C=/opt/intel_fc_80/bin/ifort
   LD=/opt/intel_fc_80/bin/ifort
   FF=-132 -O3 -c -warn noalignments
   F90F=-132 -O3 -c
   LDF=-O3 -static

endif

ifeq "$(hostname)" "panstarrs2"
   compiler = ifort
endif

ifeq "$(hostname)" "submit02"
   compiler = gfortran
endif

ifeq "$(hostname)" "mac"
   compiler = gfortran
endif
ifeq "$(hostname)" "chris"
   compiler = g95
endif

ifeq "$(hostname)" "batchelor"
	compiler =gfortran
endif

# end: special options for some hosts

# setting compiler options

ifeq "$(compiler)" "gfortran"
   FC=gfortran
   F90C=gfortran
   LD=gfortran
   FF=-ffixed-line-length-132 -O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow -Wno-align-commons -march=native -c -w -C -fdefault-real-8 -fdefault-double-8 
   F90F=-ffree-line-length-none -O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow -Wall -march=native -c -w -fdefault-real-8 -fdefault-double-8     
   LDF=-O3 -Wall -lstdc++ -finit-real=nan
   #additional new optimsation settings: -flto -fstack-arrays

  ifeq "$(warning)" "true"
     FF=-ffixed-line-length-132 -O0 -Wall  -Wno-unused -Wtabs -Wno-array-temporaries  -fno-check-array-temporaries -c -C -fdefault-real-8 -fdefault-double-8
     F90F=-ffree-line-length-none -O0 -Wall -Wtabs -Wno-array-temporaries -fno-check-array-temporaries -c -C -fdefault-real-8 -fdefault-double-8
     LDF=-O0 -Wall -lstdc++ -finit-real=nan
     DF=-fcheck=all -finit-real=snan   -Warray-bounds -Wsurprising -Wline-truncation -fbacktrace -Wunderflow  -g 
     DF90= -fcheck=all -finit-real=snan   -Warray-bounds -Wsurprising -Wline-truncation -fbacktrace -Wunderflow  -g 
     DL= 
  endif
endif

ifeq "$(compiler)" "g95"
   FC=g95
   F90C=g95
   LD=g95
   FF=-O3 -msse2 -msse -ffixed-line-length-132  -c -fbounds-check
   F90F=-O3 -msse2 -msse -ffree-line-length-huge  -c  -fbounds-check
   LDF=-O
   ifeq "$(warning)" "true"
   FF=-O0 -Wline-truncation -Wglobals -Wuninitialized -Wprecision-loss -ffixed-line-length-132  -c -fbounds-check -ftrace=full
   F90F=-O0 -Wline-truncation -Wglobals -Wuninitialized -Wprecision-loss -ffree-line-length-huge  -c  -fbounds-check -ftrace=full
  LDF=-O
  endif

endif

ifeq "$(compiler)" "ifort"
  FC=ifort
  F90C=ifort
  LD=ifort
  FF=-132 -O3 -ip -xhost -ipo  -c -vec-report0 -warn noalignments 
  F90F=-132 -O3 -ip -xhost -ipo  -c -vec-report0 
  #FF=-132 -O2 -c  -parallel 
  #F90F=-132 -O2 -c  	
  LDF=-O3  -ip -ipo
  ifeq "$(warning)" "true"
     DF= -check all -g -static -fpe0 -ftrapuv -traceback -debug all  -debug-parameters all  -warn noalignments 
     DF90= -check all -g -static -fpe0 -ftrapuv -traceback -debug all  -debug-parameters all  -warn noalignments 
     DL=-warn all -check all -g -static -fpe0 -ftrapuv -traceback -debug extended  -debug-parameters all -inline-debug-info -warn noalignments -warn nounused
  endif
endif

ifeq "$(compiler)" "ifort32"
  FC=ifort
  F90C=ifort
  LD=ifort
  FF=-132 -O3 -ip  -ipo -axWSTP  -c -vec-report0 -warn noalignments
  F90F=-132 -O3 -ip -ipo -axWSTP -c -vec-report0
  LDF=-O3 -ip -ipo -static
    ifeq "$(warning)" "true"	
       DF= -check all -g -static -fpe0 -ftrapuv -traceback -debug all  -debug-parameters all  -warn noalignments 
       DF90= -check all -g -static -fpe0 -ftrapuv -traceback -debug all  -debug-parameters all  -warn noalignments  
       DL=-warn all -check all -g -static -fpe0 -ftrapuv -traceback -debug extended  -debug-parameters all -inline-debug-info -warn noalignments -warn nounused	 
  endif
endif


OBJ = createhisto.o
all: $(OBJ) $(LIB)
	echo 'Building' ${codename}
	echo 'hostname' ${hostname}
	${LD}  ${LDF} ${DL}-o ${codename} $(MOD) $(OBJ) $(LIB)

%.o: %.f
	${FC} ${FF} ${DF}$<

%.o: %.f90 $(MOD)
	${F90C} ${F90F} ${DF90}$<


clean: clean_planete 
#       -rm *.o
#       -rm *.mod
#       cd $(DIREMPS); $(MAKE) clean

clean_planete:
	-rm *.o
	-rm *~
	-rm *.mod

