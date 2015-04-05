typedef void (*func3arg) (double, double*, double*);
typedef void (*func4arg) (double, double*, double*, double*);
typedef void (*func6arg) (double, double*, double*, double*, double*, double* );

void euler ( int neqn, double t0, double *x0, double h, func3arg f )

/*
     PURPOSE:
          perform single time step of Euler's method to approximate 
          the solution of an initial value problem (this routine works
          for a single equation as well as for a system of equations)
          

     CALLING SEQUENCE:
          euler ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp;

     xp = new double [neqn];
     
     f(t0,x0,xp);
     for ( i = 0; i < neqn; i++ )
         x0[i] += h * xp[i];
     
     delete [] xp;
}


void mod_euler ( int neqn, double t0, double *x0, double h, func3arg f )

/*
     PURPOSE:
          perform single time step of the modified Euler's method to 
          approximate the solution of an initial value problem (this 
          routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          mod_euler ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp, *xtilde;

     xp = new double [neqn];
     xtilde = new double [neqn];
     
     f(t0,x0,xp);
     for ( i = 0; i < neqn; i++ )
         xtilde[i] = x0[i] + (h/2.0) * xp[i];
     
     f(t0+h/2.0, xtilde, xp);
     for ( i = 0; i < neqn; i++ )
         x0[i] += h * xp[i];
         
     delete [] xp;
     delete [] xtilde;
}


void heun ( int neqn, double t0, double *x0, double h, func3arg f )

/*
     PURPOSE:
          perform single time step of the Heun method to approximate 
          the solution of an initial value problem (this routine works 
          for a single equation as well as for a system of equations)
          

     CALLING SEQUENCE:
          heun ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp0, *xp1, *xtilde;

     xp0 = new double [neqn];
     xp1 = new double [neqn];
     xtilde = new double [neqn];
     
     f(t0,x0,xp0);
     for ( i = 0; i < neqn; i++ )
         xtilde[i] = x0[i] + h * xp0[i];
     
     f(t0+h, xtilde, xp1);
     for ( i = 0; i < neqn; i++ )
         x0[i] += (h/2.0) * ( xp0[i] + xp1[i] );
         
     delete [] xp0;
     delete [] xp1;
     delete [] xtilde;
}


void opt_rk2 ( int neqn, double t0, double *x0, double h, func3arg f )

/*
     PURPOSE:
          perform single time step of the optimal RK2 method to approximate 
          the solution of an initial value problem (this routine works 
          for a single equation as well as for a system of equations)
          

     CALLING SEQUENCE:
          opt_rk2 ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp0, *xp1, *xtilde;

     xp0 = new double [neqn];
     xp1 = new double [neqn];
     xtilde = new double [neqn];
     
     f(t0,x0,xp0);
     for ( i = 0; i < neqn; i++ )
         xtilde[i] = x0[i] + (2.0*h/3.0) * xp0[i];
     
     f(t0+2.0*h/3.0, xtilde, xp1);
     for ( i = 0; i < neqn; i++ )
         x0[i] += ( (h/4.0)*xp0[i] + (3.0*h/4.0)*xp1[i] );
         
     delete [] xp0;
     delete [] xp1;
     delete [] xtilde;
}


void taylor_2nd ( int neqn, double t0, double *x0, double h, func4arg f )

/*
     PURPOSE:
          perform single time step of the second-order Taylor method to 
          approximate the solution of an initial value problem (this 
          routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          taylor_2nd ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of four arguments which defines the
                        right-hand side of the differential equation
                        and its first derivative with respect to the
                        independent variable;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp,
                                    double *xpt )
                           {
                              xp[0] = ...;       xpt[0] = ...;
                              xp[1] = ...;       xpt[1] = ...;
                              xp[2] = ...;       xpt[2] = ...;
                                    .                   .
                                    .                   .
                                    .                   .
                              xp[neqn-1] = ...;  xpt[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system and xpt[i] is
                        the value of the first derivative with respect to
                        the independent variable of the i-th equation in
                        the system
                        type:  func4arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp, *xpt;

     xp  = new double [neqn];
     xpt = new double [neqn];
     
     f(t0,x0,xp,xpt);
     for ( i = 0; i < neqn; i++ )
         x0[i] += ( h*xp[i] + h*h*xpt[i]/2.0 );
         
     delete [] xp;
     delete [] xpt;
}


void taylor_4th ( int neqn, double t0, double *x0, double h, func6arg f )

/*
     PURPOSE:
          perform single time step of the fourth-order Taylor method to 
          approximate the solution of an initial value problem (this 
          routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          taylor_4th ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of six arguments which defines the
                        right-hand side of the differential equation
                        and its first three derivatives with respect 
                        to the independent variable;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp,
                                    double *xpt, double *xpt2, double *xpt3 )
                           {
                              xp[0] = ...;         xpt[0] = ...;
                              xpt2[0] = ...;       xpt3[0] = ...;
                              xp[1] = ...;         xpt[1] = ...;
                              xpt2[1] = ...;       xpt3[1] = ...;
                              xp[2] = ...;         xpt[2] = ...;
                              xpt2[2] = ...;       xpt3[2] = ...;
                                    .                     .
                                    .                     .
                                    .                     .
                              xp[neqn-1] = ...;    xpt[neqn-1] = ...;
                              xpt2[neqn-1] = ...;  xpt3[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system, xpt[i] is the 
                        value of the first derivative with respect to
                        the independent variable of the right-hand side
                        of the i-th equation, xpt2[i] is the second 
                        derivative with respect to the independent variable 
                        of the right-hand side of the i-th equation and 
                        xpt3[i] is the third derivative with respect to 
                        the independent variable of the right-hand side 
                        of the i-th equation in the system
                        type:  func6arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp, *xpt, *xpt2, *xpt3;

     xp   = new double [neqn];
     xpt  = new double [neqn];
     xpt2 = new double [neqn];
     xpt3 = new double [neqn];
     
     f(t0,x0,xp,xpt,xpt2,xpt3);
     for ( i = 0; i < neqn; i++ )
         x0[i] += ( h*xp[i] + h*h*xpt[i]/2.0 + h*h*h*xpt2[i]/6.0 
                            + h*h*h*h*xpt3[i]/24.0 );
     
     delete [] xp;
     delete [] xpt;
     delete [] xpt2;
     delete [] xpt3;
}


void rk4 ( int neqn, double t0, double *x0, double h, func3arg f )

/*
     PURPOSE:
          perform single time step of the classical fourth-order Runge-Kutta 
          method to approximate the solution of an initial value problem 
          (this routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          rk4 ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{
     int i;
     double *xp, *xtilde, *k1, *k2, *k3, *k4;

     xp = new double [neqn];
     xtilde = new double [neqn];
     k1 = new double [neqn];
     k2 = new double [neqn];
     k3 = new double [neqn];
     k4 = new double [neqn];
     
     f(t0,x0,xp);
     for ( i = 0; i < neqn; i++ ) {
         k1[i] = h*xp[i];
         xtilde[i] = x0[i] + k1[i]/2.0;
     }
     
     f(t0+h/2.0, xtilde, xp);
     for ( i = 0; i < neqn; i++ ) {
         k2[i] = h*xp[i];
         xtilde[i] = x0[i] + k2[i]/2.0;
     }
     
     f(t0+h/2.0, xtilde, xp);
     for ( i = 0; i < neqn; i++ ) {
         k3[i] = h*xp[i];
         xtilde[i] = x0[i] + k3[i];
     }
     
     f(t0+h, xtilde, xp);
     for ( i = 0; i < neqn; i++ )
         k4[i] = h*xp[i];
         
     for ( i = 0; i < neqn; i++ )
         x0[i] += ( (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0 );
     
     delete [] xp;
     delete [] xtilde;
     delete [] k1;
     delete [] k2;
     delete [] k3;
     delete [] k4;

}


void ab2 ( int neqn, double t0, double *x0, double h, int *call_num,
           double *work, func3arg f )

/*
     PURPOSE:
          perform single time step of the two-step Adams-Bashforth method 
          to approximate the solution of an initial value problem (this
          routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          ab2 ( neqn, t0, x0, h, call_num, work, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          call_num      for the first call to 'ab2' and whenever the
                        value of h for the current time step is different
                        from the value of h for the previous time step,
                        set *call_num = 1
                        type:  *int
          work          array of length at least neqn; used for saving
                        information needed for subsequent calls to 'ab2'
                        type:  *double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
          call_num      after each call to 'ab2', the value of this variable
                        is incremented by 1; unless the value of h (the 
                        length of the time step) changes from one time
                        step to the next, the value of this variable should
                        not be altered
                        type:  *int
          work          contains information needed for subsequent calls
                        to 'ab2'; the contents of this array should not
                        be altered
                        type:  *double
*/

{
     int i;
     double *u, xtilde, fold, fnew;

     if ( (*call_num) == 1 ) {
/*
   if this is the first call to 'ab2', or if the value of h has been
   changed fom the value used during the previous time step
   ...   then   ...
   generate the approximation using the optimal RK2 method and save
   the first function evaluation for the next call to 'ab2'
*/

        double *xp0, *xp1, *xtilde;
        
        xp0 = new double [neqn];
        xp1 = new double [neqn];
        xtilde = new double [neqn];
     
        f(t0,x0,xp0);
        for ( i = 0; i < neqn; i++ ) {
            xtilde[i] = x0[i] + (2.0*h/3.0) * xp0[i];
            work[i] = xp0[i];
        }
     
        f(t0+2.0*h/3.0, xtilde, xp1);
        for ( i = 0; i < neqn; i++ )
            x0[i] += ( (h/4.0)*xp0[i] + (3.0*h/4.0)*xp1[i] );
         
        delete [] xp0;
        delete [] xp1;
        delete [] xtilde;
     }
     else {

/*
   use the two-step Adams-Bashforth method 
*/ 
    
        double *xp;
        xp = new double [neqn];
        
        f(t0,x0,xp);
        for ( i = 0; i < neqn; i++ ) {
            x0[i] += h*( 1.5*xp[i] - 0.5*work[i] );
            work[i] = xp[i];
        }
        
        delete [] xp;
     }
     (*call_num)++;

}


void ab4 ( int neqn, double t0, double *x0, double h, int *call_num,
           double *work, func3arg f )

/*
     PURPOSE:
          perform single time step of the four-step Adams-Bashforth method 
          to approximate the solution of an initial value problem (this
          routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          ab4 ( neqn, t0, x0, h, call_num, work, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          call_num      for the first call to 'ab4' and whenever the
                        value of h for the current time step is different
                        from the value of h for the previous time step,
                        set *call_num = 1
                        type:  *int
          work          array of length at least 3*neqn; used for saving
                        information needed for subsequent calls to 'ab4'
                        type:  *double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
          call_num      after each call to 'ab4', the value of this variable
                        is incremented by 1; unless the value of h (the 
                        length of the time step) changes from one time
                        step to the next, the value of this variable should
                        not be altered
                        type:  *int
          work          contains information needed for subsequent calls
                        to 'ab4'; the contents of this array should not
                        be altered
                        type:  *double
*/

{
     int i;
     double *u, xtilde, fold, fnew;

     if ( (*call_num) < 4 ) {
/*
   if this is one of the first three calls to 'ab4', or if the value of h 
   has been changed fom the value used during the previous time step
   ...   then   ...
   generate the approximation using the classical fourth-order Runge-Kutta
   method and save the first function evaluation for the next call to 'ab4'
*/

        double *xp, *xtilde, *k1, *k2, *k3, *k4;

        xp = new double [neqn];
        xtilde = new double [neqn];
        k1 = new double [neqn];
        k2 = new double [neqn];
        k3 = new double [neqn];
        k4 = new double [neqn];
     
        f(t0,x0,xp);
        for ( i = 0; i < neqn; i++ ) {
            k1[i] = h*xp[i];
            xtilde[i] = x0[i] + k1[i]/2.0;
            work[((*call_num)-1)*neqn+i] = xp[i];
        }
     
        f(t0+h/2.0, xtilde, xp);
        for ( i = 0; i < neqn; i++ ) {
            k2[i] = h*xp[i];
            xtilde[i] = x0[i] + k2[i]/2.0;
        }
     
        f(t0+h/2.0, xtilde, xp);
        for ( i = 0; i < neqn; i++ ) {
            k3[i] = h*xp[i];
            xtilde[i] = x0[i] + k3[i];
        }
     
        f(t0+h, xtilde, xp);
        for ( i = 0; i < neqn; i++ )
            k4[i] = h*xp[i];
         
        for ( i = 0; i < neqn; i++ )
            x0[i] += ( (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0 );
     
        delete [] xp;
        delete [] xtilde;
        delete [] k1;
        delete [] k2;
        delete [] k3;
        delete [] k4;
     }
     else {

/*
   use the four-step Adams-Bashforth method 
*/ 
    
        double *xp;
        xp = new double [neqn];
        
        f(t0,x0,xp);
        for ( i = 0; i < neqn; i++ ) {
            x0[i] += (h/24.0)*( 55.0*xp[i] - 59.0*work[2*neqn+i] 
                                + 37.0*work[neqn+i] - 9.0*work[i] );
            work[i] = work[neqn+i];
            work[neqn+i] = work[2*neqn+i];
            work[2*neqn+i] = xp[i];
        }
        
        delete [] xp;
     }
     (*call_num)++;

}


void adams_pc4 ( int neqn, double t0, double *x0, double h, int *call_num,
                 double *work, func3arg f )

/*
     PURPOSE:
          perform single time step of the Adams fourth-order 
          predictor-corrector method to approximate the solution of an 
          initial value problem (this routine works for a single equation 
          as well as for a system of equations)
          

     CALLING SEQUENCE:
          adams_pc4 ( neqn, t0, x0, h, call_num, work, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0		initial value for independent variable
                        type:  double
          x0		initial value for dependent variable
                        type:  *double
          h		length of time step
                        type:  double
          call_num      for the first call to 'adams_pc4' and whenever the
                        value of h for the current time step is different
                        from the value of h for the previous time step,
                        set *call_num = 1
                        type:  *int
          work          array of length at least 3*neqn; used for saving
                        information needed for subsequent calls to 'adams_pc4'
                        type:  *double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg


     OUTPUT:
          x0		approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
          call_num      after each call to 'adams_pc4', the value of this 
                        variable is incremented by 1; unless the value of h 
                        (the length of the time step) changes from one time
                        step to the next, the value of this variable should
                        not be altered
                        type:  *int
          work          contains information needed for subsequent calls
                        to 'adams_pc4'; the contents of this array should not
                        be altered
                        type:  *double
*/

{
     int i;
     double *u, xtilde, fold, fnew;

     if ( (*call_num) < 4 ) {
/*
   if this is one of the first three calls to 'adams_pc4', or if the value 
   of h has been changed fom the value used during the previous time step
   ...   then   ...
   generate the approximation using the classical fourth-order Runge-Kutta
   method and save the first function evaluation for the next call to 
   'adams_pc4'
*/

        double *xp, *xtilde, *k1, *k2, *k3, *k4;

        xp = new double [neqn];
        xtilde = new double [neqn];
        k1 = new double [neqn];
        k2 = new double [neqn];
        k3 = new double [neqn];
        k4 = new double [neqn];
     
        f(t0,x0,xp);
        for ( i = 0; i < neqn; i++ ) {
            k1[i] = h*xp[i];
            xtilde[i] = x0[i] + k1[i]/2.0;
            work[((*call_num)-1)*neqn+i] = xp[i];
        }
     
        f(t0+h/2.0, xtilde, xp);
        for ( i = 0; i < neqn; i++ ) {
            k2[i] = h*xp[i];
            xtilde[i] = x0[i] + k2[i]/2.0;
        }
     
        f(t0+h/2.0, xtilde, xp);
        for ( i = 0; i < neqn; i++ ) {
            k3[i] = h*xp[i];
            xtilde[i] = x0[i] + k3[i];
        }
     
        f(t0+h, xtilde, xp);
        for ( i = 0; i < neqn; i++ )
            k4[i] = h*xp[i];
         
        for ( i = 0; i < neqn; i++ )
            x0[i] += ( (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0 );
     
        delete [] xp;
        delete [] xtilde;
        delete [] k1;
        delete [] k2;
        delete [] k3;
        delete [] k4;
     }
     else {

/*
   use the Adams fourth-order predictor-corrector method 
*/ 
    
        double *xp0, *xtilde, *xp1;
        xp0 = new double [neqn];
        xtilde = new double [neqn];
        xp1 = new double [neqn];
        
        f(t0,x0,xp0);
        for ( i = 0; i < neqn; i++ ) 
            xtilde[i] = x0[i] + (h/24.0)*( 55.0*xp0[i] - 59.0*work[2*neqn+i] 
                                           + 37.0*work[neqn+i] - 9.0*work[i] );
        
        f(t0+h,xtilde,xp1);
        for ( i = 0; i < neqn; i++ ) {
            x0[i] += (h/24.0)*( 9.0*xp1[i] + 19.0*xp0[i] - 5.0*work[2*neqn+i] 
                                + work[neqn+i] );
            work[i] = work[neqn+i];
            work[neqn+i] = work[2*neqn+i];
            work[2*neqn+i] = xp0[i];
        }
        
        delete [] xp0;
        delete [] xtilde;
        delete [] xp1;
     }
     (*call_num)++;

}
