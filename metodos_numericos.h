/**
 * @file metodos_numericos.h
 * @author Tanil Izquierdo Cordova.
 * @brief Libreria para el calculo de metodos numericos iterativos.
 * @version 1
 * @date 2022-10-31
 * 
 * @copyright GPL V3
 * 
 */

#ifndef _FERMI_INCLUDED_
#define _FERMI_INCLUDED_
#include "atmsp.h"

// Valor utilizado por estimadoras Texas Instruments para derivar por medio del metodo del cociente de diferencia simetrica.
#define H 0.001

// Numero de divisiones de la integral
#define PASO_INTEGRAL 10000000

/**
 * @brief Funcion que estima la derivada de otra funcion con el metodo del cociente de diferencia simetrica.
 * 
 * @param f Funcion a derivar expresada en C.
 * @param x punto de la funcion que se desea derivar.
 * @param orden Orden de la derivada.
 * @param h Cantidad pequeña necesaria para poder realizar la operacion de cociente de diferencia simetrica.
 * @return Aproximacion de la derivada dada f(x): (f^[orden](x)). 
 */
double derivar(double (* f)(double), double x, int orden, double h){
    if (orden > 1) {
        return (derivar(f, x+h, orden-1, h) - derivar(f, x-h, orden-1, h)) / (2.0*h);
    }
    return (f(x+h) - f(x-h)) / (2.0*h);
}

/**
 * @brief Funcion que estima la derivada de otra funcion con el metodo del cociente de diferencia simetrica.
 * 
 * @param expresion Expresion matematica representativa de la funcion. f(x) = expresion.
 * @param x punto de la funcion que se desea derivar.
 * @param orden Orden de la derivada.
 * @param h Cantidad pequeña necesaria para poder realizar la operacion de cociente de diferencia simetrica.
 * @return Aproximacion de la derivada dada f(x): (f^[orden](x)). 
 */
double derivar(std::string expresion, double x, int orden, double h) {
    ATMSB <double> f;
    ATMSP <double> parser;
    parser.parse(f, expresion, "x");

    if (orden > 1) {
        return (derivar(expresion, x+h, orden-1, h) - derivar(expresion, x-h, orden-1, h)) / (2.0*h);
    }

    f.var[0] = x+h;
    double fpx = f.run();

    f.var[0] = x-h;
    double fmx = f.run();

    return (fpx - fmx) / (2.0*h);
}

/**
 * @brief Funcion que estima la derivada de una funcion que recibe dos argumetos con respecto de una variable.
 * 
 * @param f Funcion a derivar expresada en C.
 * @param x Valor de x en el que se desea derivar.
 * @param y Valor en y en el que se desea derivar.
 * @param variable La variable respecto a la cual se desea derivar.
 * @param h antidad pequeña necesaria para poder realizar la operacion de cociente de diferencia simetrica.
 * @return  Aproximacion de la derivada dada f(x, y): (f^(x, y)). 
 */
double derivar(double (* f)(double, double), double x, double y, char variable, double h){
    double resultado = 0;
    if (variable == 'x') {
        resultado = (f(x+h, y) - f(x-h, y)) / (2.0*h);
    }
    else if (variable == 'y') {
        resultado = (f(x, y+h) - f(x, y-h)) / (2.0*h);
    }
    return resultado;
}

/**
 * @brief Funcion que estima la derivada de una funcion que recibe dos argumetos con respecto de una variable.
 * 
 * @param expresion Expresion matematica representativa de la funcion. f(x, y) = expresion.
 * @param x Valor de x en el que se desea derivar.
 * @param y Valor en y en el que se desea derivar.
 * @param variable La variable respecto a la cual se desea derivar.
 * @param h Cantidad pequeña necesaria para poder realizar la operacion de cociente de diferencia simetrica.
 * @return Aproximacion de la derivada dada f(x, y): (f^(x, y)). 
 */
double derivar(std::string expresion, double x, double y, char variable, double h){
    double resultado = 0;
    ATMSB <double> f;
    ATMSP <double> parser;
    parser.parse(f, expresion, "x, y");
    if (variable == 'x') {
        f.var[0] = x+h;
        f.var[1] = y;
        double fpx = f.run();

        f.var[0] = x-h;
        double fmx = f.run();

        resultado = (fpx - fmx) / (2.0*h);
    }
    else if (variable == 'y') {
        f.var[0] = x;
        f.var[1] = y+h;
        double fpy = f.run();

        f.var[1] = y-h;
        double fmy = f.run();

        resultado = (fpy - fmy) / (2.0*h);
    }
    return resultado;
}

double error_integracion_simpson(double (*f)(double), double a, double b) {
    double xi = a + ((b-a)/2);
    return -(pow(b-a, 5)/6480) * derivar(f, xi, 4, H);
}

/**
 * Funcion que integra en base a la regla trapezoidal aproximando la integral.
 * La funcion sobrecargada `integrar(double f(double) | string, a, b)` es superior y deberia de ser utilizada en lugar de esta.
 * @param f Funcion a integrar definida en C.
 * @param a Limite inferior.
 * @param b Limite superior.
 * @param n Numero de subdivisiones para aproximar la integral, se incluye la constante `PASO_INTEGRAL` para llegar a una precision aceptable.
 * 
 * @return Aproximacion de la integral.
 */
double integrar(double (* f)(double), double a, double b, double n) {
    double sigma = 0;
    double h = (b-a)/n;

    for (int k = 1; k < n-1; k++) {
        sigma += f(a + k * h);
    }

    return h * ((f(a)/2) + sigma + (f(b)/2));
}

/**
 * Función que integra en base a la regla de 3/8 de Simpson, una variación de las formulas de Newton-Cotes.
 * https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule
 * @param f Funcion a integrar definida en C.
 * @param a Limite inferior.
 * @param b Limite superior. 
 * 
 * @return Aproximación de la integral.
 */
double integrar(double (* f)(double), double a, double b){
    //return ((b-a)/6) * (f(a) +4*f((a+b)/2) + f(b));
    return (((b-a)/8) * (f(a) + 3*f((2*a+b)/3) + 3*f((a+2*b)/3) + f(b)));
    //+ error_integracion_simpson(f, a, b);
}

/**
 * @brief Funcion que integra en base a la regla de 3/8 de Simpson, una variación de las formulas de Newton-Cotes.
 * https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule.
 * 
 * @param expresion Expresion representativa de la funcion. f(x) = expresion.
 * @param a Limite inferior.
 * @param b Limite superior.
 * 
 * @return Aproximación de la integral.
 */
double integrar(std::string expresion, double a, double b){
    ATMSB <double> f;
    ATMSP <double> parser;
    parser.parse(f, expresion, "x");

    f.var[0] = a;
    double fa = f.run();

    f.var[0] = b;
    double fb = f.run();

    f.var[0] = (2*a+b)/3;
    double famul = f.run();

    f.var[0] = (a+2*b)/3;
    double fbmul = f.run();

    return ((b-a)/8) * (fa + 3*famul + 3*fbmul + fb);
}

double integracion_trapecio(double (*f)(double), double a, double b, int n) {
    double paso = (b-a)/n;
    double suma = 0;
    for (double xi = a+paso; xi < b; xi += paso) {
        suma += f(xi);
    }

    return (b-a)*(f(a) + 2*suma + f(b))/(2*n);
}

/**
 * @brief Funcion para estimar raices reales mediante el metodo de Newton-Raphson.
 * 
 * @param u Funcion definida en C.
 * @param v Funcion definida en C.
 * @param x Valor inicial (x_0).
 * @param y Valor inicial (y_0).
 * @param it Iteraciones deseadas.
 * @return Arreglo de longitud 2 con los valores de "x" en la posicion 0 y de "y" en la posicion 1 para la estimacion de las raices reales.
 */
double * newton_raphson_nolineal(double (*u)(double, double), double (*v)(double, double), double x, double y, int it) {
    for (int i=0; i < it; i++) {
        x = x - (u(x, y)*derivar(v, x, y, 'y', H) - v(x, y)*derivar(u, x, y, 'y', H))
                /((derivar(u, x, y, 'x', H)*derivar(v, x, y, 'y', H)) - (derivar(u, x, y, 'y', H)*derivar(v, x, y, 'x', H)));

        y = y - (v(x, y)*derivar(u, x, y, 'x', H) - u(x, y)*derivar(v, x, y, 'x', H))
                /((derivar(u, x, y, 'x', H)*derivar(v, x, y, 'y', H)) - (derivar(u, x, y, 'y', H)*derivar(v, x, y, 'x', H)));
    }

    static double resultados[2];
    resultados[0] = x;
    resultados[1] = y;
    return resultados;
}

/**
 * @brief Funcion para estimar raices reales mediante el metodo de Newton-Raphson.
 * 
 * @param expresion_u Expresion de representativa de la funcion u. u(x, y) = expresion.
 * @param expresion_v Expresion de representativa de la funcion v. v(x, y) = expresion.
 * @param x Valor inicial (x_0).
 * @param y Valor inicial (y_0).
 * @param it Iteraciones deseadas.
 * @return Arreglo de longitud 2 con los valores de "x" en la posicion 0 y de "y" en la posicion 1 para la estimacion de las raices reales.
 */
double * newton_raphson_nolineal(std::string expresion_u, std::string expresion_v, double x, double y, int it) {
    ATMSB <double> u;
    ATMSB <double> v;
    ATMSP <double> parser;
    parser.parse(u, expresion_u, "x, y");
    parser.parse(v, expresion_v, "x, y");
    double uxy;
    double vxy;
    
    for (int i=0; i < it; i++) {
        u.var[0] = x;
        u.var[1] = y;
        uxy = u.run();

        v.var[0] = x;
        v.var[1] = y;
        vxy = v.run();

        x = x - (uxy*derivar(expresion_v, x, y, 'y', H) - vxy*derivar(expresion_u, x, y, 'y', H))
                /((derivar(expresion_u, x, y, 'x', H)*derivar(expresion_v, x, y, 'y', H)) - (derivar(expresion_u, x, y, 'y', H)*derivar(expresion_v, x, y, 'x', H)));

        u.var[0] = x;
        uxy = u.run();

        v.var[0] = x;
        vxy = v.run();

        y = y - (vxy*derivar(expresion_u, x, y, 'x', H) - uxy*derivar(expresion_v, x, y, 'x', H))
                /((derivar(expresion_u, x, y, 'x', H)*derivar(expresion_v, x, y, 'y', H)) - (derivar(expresion_u, x, y, 'y', H)*derivar(expresion_v, x, y, 'x', H)));
    }

    static double resultados[2];
    resultados[0] = x;
    resultados[1] = y;
    return resultados;
}

#endif