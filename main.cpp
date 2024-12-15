#include <iostream>
#include <cmath>
#include <boost/program_options.hpp>
#include <tabulate/tabulate.hpp>
#include "fermi.h"

namespace po = boost::program_options;

int main(int narg, char ** args) {
    po::options_description descripciones ("Tanil Izquierdo Cordova.\nPrograma para el calculo de métodos numericos iterativos");
    descripciones.add_options() ("help,h", "Muestra este mensaje de ayuda.")
                                ("ayuda", "Muestra este mensaje de ayuda.")
                                ("metodo,m", po::value<std::string>()->default_value(""), "Metodo numerico que se desea utilizar.");
    po::variables_map argumentos;
    po::store(po::parse_command_line(narg, args, descripciones), argumentos);
    po::notify(argumentos);

    if (argumentos.count("help") || argumentos.count("ayuda")) {
        std::cout << descripciones << "\n";
    }

    if (argumentos.count("metodo")) {
        std::cout << argumentos["metodo"].as<std::string>() << "\n";
    }

    //Table nr_nolineal;

    //nr_nolineal.format().multi_byte_characters(true);

    //nr_nolineal[0].format()
    //.padding_top(1)
    //.padding_bottom(1)
    //.font_align(FontAlign::center)
    //.font_style({FontStyle::underline})
    //.font_background_color(Color::green);

    //std::cout << nr_nolineal << std::endl;

    if (__cplusplus == 202101L) std::cout << "C++23";
    else if (__cplusplus == 202002L) std::cout << "C++20";
    else if (__cplusplus == 201703L) std::cout << "C++17";
    else if (__cplusplus == 201402L) std::cout << "C++14";
    else if (__cplusplus == 201103L) std::cout << "C++11";
    else if (__cplusplus == 199711L) std::cout << "C++98";
    else std::cout << "pre-standard C++." << __cplusplus;
    std::cout << "\n";


    std::cout << "fokin did it" << "\n";

    return EXIT_SUCCESS;
}
