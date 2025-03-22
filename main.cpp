#include "Elem4.h"


Elem4::Elem4() {
}

int main() {
    Elem4 mes;
    
    //zakomentuj pozosta³e aby wybraæ

    mes.loadData("Test1_4_4.txt");
    //mes.loadData("Test2_4_4_MixGrid.txt");
    //mes.loadData("Test3_31_31_kwadrat.txt");


    mes.printData();
    //wybor - wybranie ilosci punktow ca³kowania
    //zmiana 2/3/4
    int wybor = 2;


    mes.Pochodne(wybor);
    mes.Obliczanie(wybor);
    mes.MacierzC(wybor);

    mes.PrintGlobalH();
    mes.PrintWektorP();

    mes.Symulacja();

    return 0;
}