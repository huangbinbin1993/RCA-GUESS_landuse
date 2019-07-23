///////////////////////////////////////////////////////////////////////////////////////
// FLINK.H
// Header file for C++ source files to be linked to a Fortran executable
// Author: Ben Smith, Lund University / Benjamin.Smith@nateko.lu.se
// Version date: 2004-08-10

///////////////////////////////////////////////////////////////////////////////////////
// SYNOPSIS
//
// This header file provides preprocessor macros that help to make C++ functions
// available as subroutines to a fortran program.
//
// NOTE: NO FUNCTIONALITY is provided to make fortran subroutines available to C++
//
// 1. To declare (and define) a C++ function to be called as a subroutine from Fortran:
//
//    FROUTINE( name )( args ) {
//        // body of function
//    }
//
//    e.g.  In C++ code:      FROUTINE(greeting)() {
//                                printf("Hello world!");
//                            }
//            
//          In Fortran code:  call greeting()
//
//    IMPORTANT: alphabetic characters in 'name' must be all lower-case
//
// 2. To declare an argument with the Fortran INTEGER type:
//
//    FINTEGER name
//
//    e.g.  In C++ code:      FROUTINE(assign)(FINTEGER to,FINTEGER from) {
//                                 FGET(to)=FGET(from)
//                            }
//
//          In Fortran code:  integer a,b
//                            b=12
//                            call assign(a,b) !sets a to 12
//
// 3. To declare an argument with the Fortran REAL type:
//
//    FREAL name
//
// 4. To declare an argument with the Fortran DOUBLE PRECISION type:
//
//    FDOUBLE name
//
// 5. To declare array arguments of the above types:
//
//    FINTEGERARR name
//    FREALARR name
//    FDOUBLEARR name
//
//    e.g.  In C++ code:      FROUTINE(init)(FINTEGERARR a,FINTEGER n) {
//                                 int i;
//                                 for (i=1;i<=n;i++)
//                                      FGETARR(a,i)=0;
//                            }
//
//          In Fortran code:  integer a(10)
//                            call init(a,10) !fills array a with zeroes
//
// 6. To reference a fortran argument as part of a C++ expression (whether for
//    assignment or access):
//
//    FGET( name )
//
//    e.g.  In C++ code:      FROUTINE(getarea)(FREAL area,FREAL radius) {
//                                 FGET(area)=3.142*FGET(from)*FGET(from);
//                            }
//
//          In Fortran code:  real a,r
//                            radius=2.5
//
//                            !(set a to area of a circle with radius 2.5)
//                            call getarea(a,r) 
//
// 7. To reference elements of a fortran array argument (whether for assignment or
//    access):
//
//    FGETARR( name , index )
//
//    NOTE THAT array indices are 1-based (as in Fortran), not 0-based (as normally in
//    C/C++)
//
//    e.g.  In C++ code:      FROUTINE(getareas)(FREALARR area,FREALARR radius,
//                                               FINTEGER nelement) {
//                                 int i;
//                                 for (i=1;i<=FGET(nelement);i++)
//                                      FGETARR(area,i)=3.142*FGETARR(radius,i)*
//                                                      FGETARR(radius,i);
//                            }
//
//          In Fortran code:  real a(10),r(10)
//                            call getradii(r) !get radii for 10 circles from somewhere
//                            call getareas(a,r,10)


///////////////////////////////////////////////////////////////////////////////////////
// HOW TO BUILD AN EXECUTABLE FROM FORTRAN AND C++ SOURCE CODE FILES
//
// 1. The "main" should be in Fortran (otherwise you would need a way to make fortran
//    subroutines available for calling from C++, and this file does not help with
//    that).
//
// 2. There may be one or more source code files in C++, these must not contain a
//    'main' (if you want an existing C++ main function to be the de facto main
//    program, rename it to something else, declared with an FROUTINE, and call it from
//    your Fortran program).
//
// 3. Any C++ files containing functions to be available for calling from Fortran
//    should #include this header file
//
// 4. First compile the C++ file(s) to object files, e.g.
//
//    CC -fast -c mycode.cpp
//
//    (this builds the object file mycode.o)
//
// 5. Compile the fortran program and link it to the C++ object files, e.g.
//
//    f77 -o myexe -fast mycode.o mymain.f
//
//    (builds the executable, myexe)


#define FROUTINE(name) extern "C" void name##_
#define FINTEGER int*
#ifdef SINGLE_PRECISION
#define FREAL float*
#else
#define FREAL double*
#endif
#define FDOUBLE double*
#define FINTEGERARR int*
#define FREALARR float*
#define FDOUBLEARR double*
#define FGET(value) *value
#define FGETARR(name,index) name[index-1]
