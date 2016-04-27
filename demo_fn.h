     struct quadratic_params
       {
         double a, b, c;
       };
     
     double quadratic2 (double x, void *params);
     double quadratic_deriv (double x, void *params);
     void quadratic_fdf (double x, void *params, 
                         double *y, double *dy);



     struct qubic_params
       {
         double a, b, c, d;
       };
     
     double qubic (double x, void *params);
     double qubic_deriv (double x, void *params);
     void qubic_fdf (double x, void *params, 
                    double *y, double *dy);

     struct quartic_params
       {
         double a, b, c, d, e;
       };
     
     double quartic (double x, void *params);
     double quartic_deriv (double x, void *params);
     void quartic_fdf (double x, void *params, 
                    double *y, double *dy);
