#include <iostream>
#include<vector>
#include <cmath>
#include "CMyVektor.cpp"

class CMyMatrix {
private:

    unsigned int spalten;
    unsigned int zeilen;
    std::vector<std::vector<double>> matrix;

public:
    CMyMatrix(unsigned int zeilen, unsigned int spalten) : spalten(spalten), zeilen(zeilen), matrix(zeilen, std::vector<double>(spalten, 0.0)){};

    void setMatrixComponent(unsigned int zeile, unsigned int spalte, double wert) {//double statt int!!!!
        matrix[zeile][spalte] = wert;
    }

    double getMatrixComponent(unsigned int  zeile, unsigned int  spalte) const {
        return matrix[zeile][spalte];
    }


    CMyMatrix invers(){
        CMyMatrix *mtx = this;
        if(mtx->zeilen != 2 || mtx->spalten != 2) {std::cout<<"verschiedene dimensionen Inverse kann nicht erzeugt werden"<<std::endl;}
        double a = mtx->getMatrixComponent(0, 0);
        double b = mtx->getMatrixComponent(0, 1);
        double c = mtx->getMatrixComponent(1, 0);
        double d = mtx->getMatrixComponent(1, 1);

        //det prüfung A!=0
        double det = (a*d)-(b*c);
        if(det == 0)
        {
            CMyMatrix n(0, 0);
            return n;
        }

        CMyMatrix inverse(2,2);

        inverse.setMatrixComponent(0, 0, d/det);
        inverse.setMatrixComponent(0, 1, -b/det);
        inverse.setMatrixComponent(1, 0, -c/det);
        inverse.setMatrixComponent(1, 1, a/det);

        return inverse;
    }

    int getDimension() {
        return this->zeilen;
    }

    int getSpalten() {
        return this->spalten;
    }




    CMyVektor operator*(const CMyVektor& vektor) {
        //debugging ausgaben
        if (this->spalten != vektor.getVektorDimension()) {
            std::cout<<"vektor dim: "<<vektor.getVektorDimension()<<std::endl;
            std::cout<<"matrix spalten: "<<this->spalten<<std::endl;
            std::cout<<"matrix zeilen: "<<this->zeilen<<std::endl;

           std::cout<<"Matrix und Vektor sind nicht kompatibel für Multiplikation."<<std::endl;
        }

        CMyVektor ergebnisVektor(this->zeilen);  //vektor dimension == matrix dimension
        for (int i = 0; i < this->zeilen;i++) {
            double summe = 0.0;
            for (int j = 0; j < this->spalten; j++) {
                summe += this->matrix[i][j] * vektor.getVektorComponent(j); //mtx zeile mal vektorkomponente
            }
            ergebnisVektor.setVektorComponent(i, summe);
        }

        return ergebnisVektor;
    }



    void print() {
        std::cout<<"matrix spalten: "<<this->spalten<<std::endl;
        std::cout<<"matrix zeilen: "<<this->zeilen<<std::endl;

    }



    //----------------------------------------------------------------------------------------------------------------//

    CMyMatrix jacobi(CMyVektor x, CMyVektor (*funktion)(CMyVektor x)) {
        double h = 1e-4; //10^-4

        CMyVektor fx = funktion(x);
        CMyMatrix jacobi =CMyMatrix(fx.getVektorDimension(), x.getVektorDimension());

        for(int i=0; i<x.getVektorDimension(); i++) {
            CMyVektor xplusH = x;
            xplusH.setVektorComponent(i,x.getVektorComponent(i)+h);
            CMyVektor fxH = funktion(xplusH);
            for(int j=0; j<fx.getVektorDimension(); j++) {

                double diff= fxH.getVektorComponent(j)-fx.getVektorComponent(j); //differenz der vektoren bilden
                jacobi.setMatrixComponent(j,i,diff/h);
            }
        }
        return jacobi;
    }

    //----------------------------------------------------------------------------------------------------------------//


    CMyVektor newtonVerfahren(CMyVektor stelle, CMyVektor (*funktion)(CMyVektor)) {
        const double genauigkeit = 1e-5;
        CMyVektor x = stelle;
        int schritte = 0;

        do {
            CMyVektor fx = funktion(x);  // f(x)
            double betrag = fx.getVektorLength();

            if (betrag < genauigkeit) {
                std::cout<<"betrag erreicht"<<std::endl;
                std::cout<<"betrag: "<<fx.getVektorLength()<<std::endl;

                break;
            }

            CMyMatrix jacobiMtx = jacobi(x, funktion); // Jacobimatrix
            CMyMatrix jacobiMtxInv = jacobiMtx.invers();       // Inverse der Jacobimatrix
            jacobiMtxInv.print();

            CMyVektor delta = jacobiMtxInv * (-1 * fx);      // delta = -J^(-1) * f(x)
            x = x + delta;                            // x = x + delta

            std::cout<<"schritt: "<<schritte<<std::endl;
            std::cout<<"betrag: "<<fx.getVektorLength()<<std::endl;
            schritte++;
        }while (schritte < 50);

        return x;
    }

};





