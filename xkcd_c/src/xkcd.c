//
// Created by Onionslime on 3/15/2022.
//

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include "nmath.h"

SEXP rxkcd(SEXP n, SEXP sd)
{
    if (! isInteger(n))
        error("'n' must be type integer");
    if (! isReal(sd))
        error("'sd' must be type double");
    if (sd <= 0)
        error("'sd' must be positive");
    int c_n = INTEGER(n)[0];
    double mean = 0;
    double c_sd = REAL(sd)[0];

    SEXP result;
    PROTECT(result = allocVector(REALSXP, c_n));



    for (int i = 0; i < c_n; i++) {
        double temp = rnorm( mean, c_sd);
        double temp2 = dnorm(temp, mean, c_sd, 0);
        REAL(result)[i] = runif(mean, temp2);
    }


    UNPROTECT(1);
    return result;
}



double dxkcd_unit(double x, double sd, int log_c, int SEP){
    // Constants
    double pi = M_PI;
    double sq2psd = 0;
    double up_bd = 0;
    //double log2sd2 = 0;
    //double log2psd2 = 0;


    double y = 0;
    sq2psd = sqrt(2 * pi) * sd;

    up_bd = 1 / sq2psd;

    //log2sd2 = log(2 * sd * sd);

    //log2psd2 = log2sd2 + log(pi);




    /*if(x > 0 && x < up_bd ) {

        if (log_c == 0 && SEP == 0) {
            y = 2 * sqrt(((-1 / 2) * log2psd2 - log(x)) * 2 * sd * sd);

        } else if (log_c == 0 && SEP == 1) {
            y = 2 * sqrt((-2) * sd * sd * log1p(-sq2psd * x));
        } else if (log_c == 1 && SEP == 0) {
            y = log(2) + 1 / 2 * (log((-1 / 2) * log2sd2 - 1 / 2 * log(pi) - log(x)) + log2sd2);
        } else if (log_c == 1 && SEP == 1) {
            y = log(2) + 1 / 2 * (log(-log1p(-sq2psd * x)) + log2sd2);
        }

    }*/


    if(x >0 && x < up_bd){
        double const1=-log(sqrt(2*pi)*sd*x);
        double const2=log(2*sd*sd);
        if (log_c==0 && SEP == 0){
            y= 2*sqrt(2)*sd*sqrt(-log(sqrt(2*pi)*x*sd));
        }

        else if(log_c==1 && SEP == 0){
            y= log(2)+0.5*const2+0.5*log(const1);
        }

        else if(log_c==0 && SEP == 1){
            y= 2*sqrt(2)*sd*sqrt(-log1p(-sqrt(2*pi)*sd*x));
        }

        else {
            y= log(2) + 0.5*const2+0.5*log(-log1p(-sqrt(2*pi)*sd*x));
        }
    }

    if (x == 0) {
        if (log_c == 0 && SEP == 0) { y = INFINITY; }
        else if (log_c == 0 && SEP == 1) { y = 0; }
        else if (log_c == 1 && SEP == 0) { y = INFINITY; }
        else if (log_c == 1 && SEP == 1) { y = -INFINITY; }
    }

    if(x == up_bd) {
        if (log_c == 0 && SEP == 0) { y = 0; }
        else if (log_c == 0 && SEP == 1) { y = INFINITY; }
        else if (log_c == 1 && SEP == 0) { y = -INFINITY; }
        else if (log_c == 1 && SEP == 1) { y = -INFINITY; }
    }

    return y;

}


SEXP dxkcd(SEXP r_x, SEXP r_sd, SEXP r_log, SEXP r_swap_end_point, SEXP r_n, SEXP r_sd_n) {
    // transfer type
    double *x = REAL(r_x);
    double *sd = REAL(r_sd);
    int log_c = INTEGER(r_log)[0];
    int SEP = INTEGER(r_swap_end_point)[0];
    int n = INTEGER(r_n)[0];
    int sd_n = INTEGER(r_sd_n)[0];


    SEXP output;
    PROTECT(output = allocVector(REALSXP, n));

    if(sd_n == 1){
        for (int i = 0; i < n; i++){
            REAL(output)[i] = dxkcd_unit(x[i], sd[0], log_c, SEP);

        }

    } else {

        for (int i = 0; i < n; i++){
            REAL(output)[i] = dxkcd_unit(x[i], sd[i], log_c, SEP);
        }

    }

    UNPROTECT(1);
    return output;


}



double pxkcd_unit(double q, double sd, int log_p, int SEP) {
    // Constants
    double pi = M_PI;

    //printf("%s", "ip1");
    double up_bd = 1/(sqrt(2*pi)*sd);
    double dy=dxkcd_unit(q,sd,0,0)/2;
    double dw=dxkcd_unit(q,sd,0,1)/2;
    double temp = 2*sqrt(2*sd*sd*sd*sqrt(2*pi))*2/3;

    double result = 0;
    printf("%s \n", "0 return");
    if (isnan(q)||isnan(sd)||isinf(sd)||sd<=0){
        result = NAN;
        error("cannot have NAN input");
        return result;
    }

    if (q>=up_bd && log_p==0){
        printf("%s \n", "first return");
        result =  1;
        return result;
    }

    if (q<=0 && log_p ==0){
        printf("%s \n", "second return");
        result =  0;
        return result;
    }

    if (q>=up_bd && log_p ==1){
        printf("%s \n", "3 return");
        result =  0;
        return result;
    }

    if (q<=0 && log_p ==1){
        printf("%s \n", "4 return");
        result = -INFINITY;
        return result;
    }



    if (log_p==1 && SEP == 0){

        result=log(2)+pnorm(-dy,0,sd,1,1)+log1p( q*dy/pnorm(-dy,0,sd,1,0) );


    }

    if(log_p==0 && SEP == 0){

        result=2*pnorm(-dy,0,sd,1,0)+2*q*dy;

    }

    if(log_p==1 && SEP == 1 && q<= pow(10,-10)){

        result=(log(temp)+3/2*log(q));

    }

    if(log_p==1 && SEP == 1  && q>pow(10,-10)){

        result = log(1-2*pnorm(-dw,0,sd,1,0)-2*(up_bd-q)*dw);

    }

    if(log_p==0 && SEP == 1 &&  q<= pow(10,-10)){

        result=temp*pow(q,1.5);


    }

    if(log_p==0 && SEP == 1 && q> pow(10,-10)){

        result= 1-(2*pnorm(-dw,0,sd,1,0)+2*(up_bd-q)*dw);

    }
    return result;
}









SEXP pxkcd(SEXP r_x, SEXP r_sd, SEXP r_log, SEXP r_swap_end_point, SEXP r_n, SEXP r_sd_n) {
    // transfer type
    double *x = REAL(r_x);
    double *sd = REAL(r_sd);
    int log_c = INTEGER(r_log)[0];
    int SEP = INTEGER(r_swap_end_point)[0];
    int n = INTEGER(r_n)[0];
    int sd_n = INTEGER(r_sd_n)[0];

    SEXP output;
    PROTECT(output = allocVector(REALSXP, n));

    if(sd_n == 1){

        for (int i = 0; i < n; i++){

            REAL(output)[i] = pxkcd_unit(x[i], sd[0], log_c, SEP);

        }

    } else {



        for (int i = 0; i < n; i++){

            REAL(output)[i] = pxkcd_unit(x[i], sd[i], log_c, SEP);


        }

    }

    UNPROTECT(1);
    return output;


}










SEXP qxkcd(SEXP r_p, SEXP r_sd, SEXP r_log, SEXP r_swap_end_point, SEXP r_n, SEXP r_sd_n) {
    // transfer type
    double *p = REAL(r_p);
    double *sd = REAL(r_sd);
    int log_c = INTEGER(r_log)[0];
    int SEP = INTEGER(r_swap_end_point)[0];
    int n = INTEGER(r_n)[0];
    int sd_n = INTEGER(r_sd_n)[0];

    double pi = M_PI;


    double func(double q, double p, double sd, int log_p, int SEP)
    {
        return pxkcd_unit(q, sd, log_p, SEP)- p;
    }


    double bisection(double low,double up, double p, double sd, int log_p, int SEP)
    {
        double e=0.001;
        double c = 0;

        //if(func(low, p ,sd, log_p, SEP) * func(up, p ,sd, log_p, SEP) >= 0)
        //{
        //    printf("%f \n",func(low, p ,sd, log_p, SEP));
        //    printf("%f \n", func(up, p ,sd, log_p, SEP));
        //    error("Incorrect a and b");
        //}

        c = low;

        while ((up-low) >= e)
        {
            c = (low+up)/2;
            if (func(c, p ,sd, log_p, SEP) == 0.0){
                printf("Root = %lf\n",c);
                break;
            }
            else if (func(c, p ,sd, log_p, SEP)*func(low, p ,sd, log_p, SEP) < 0){
                //printf("Root = %lf\n",c);
                up = c;
            }
            else{
                //printf("Root = %lf\n",c);
                low = c;
            }
        }
        return c;
    }

    //printf("%s \n", "p1");

    SEXP output;
    PROTECT(output = allocVector(REALSXP, n));
    if(sd_n == 1) {
        for (int i = 0; i < n; ++i) {
            double low = 0;
            double up = 1 / (sqrt(2 * pi) * sd[0]);
            REAL(output)[i] = bisection(low, up, p[i], sd[0], log_c, SEP);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            double low = 0;
            double up = 1 / (sqrt(2 * pi) * sd[i]);
            REAL(output)[i] = bisection(low, up, p[i], sd[i], log_c, SEP);
        }
    }

    UNPROTECT(1);
    return output;


}