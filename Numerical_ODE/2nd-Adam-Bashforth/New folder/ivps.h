typedef void (*func3arg) (double, double*, double*);
typedef void (*func4arg) (double, double*, double*, double*);
typedef void (*func6arg) (double, double*, double*, double*, double*, double* );

void euler ( int, double, double*, double, func3arg );
void mod_euler ( int, double, double*, double, func3arg );
void heun ( int, double, double*, double, func3arg );
void opt_rk2 ( int, double, double*, double, func3arg );
void taylor_2nd ( int, double, double*, double, func4arg );
void taylor_4th ( int, double, double*, double, func6arg );
void rk4 ( int, double, double*, double, func3arg );
void ab2 ( int, double, double*, double, int*, double*, func3arg  );
void ab4 ( int, double, double*, double, int*, double*, func3arg  );
void adams_pc4 ( int, double, double*, double, int*, double*, func3arg  );