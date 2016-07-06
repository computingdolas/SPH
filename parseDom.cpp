#include "Parser.h"
#include "VTKWriter.h"

int main(){
    Parser p("bubbleWithTank");
    p.readParameters();
    p.readInputConfiguration();

    VTKWriter v("tester");
    v.write(mass,position,velocity,);

}
