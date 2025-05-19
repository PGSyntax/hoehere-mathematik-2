//
// Created by legue on 03.05.2025.
//
#include "CMyMatrix.cpp"
#include <iomanip>


class CDGLSolver {
private:
    // Funktionspointer für DGL-System
    CMyVektor (*f_DGL_System)(CMyVektor y, double x);
    // Funktionspointer für DGL höherer Ordnung
    double (*f_DGL_nterOrdnung)(CMyVektor y, double x);
    bool isHoehere;

    int schritteHEun = 0;


    CMyVektor ableitungen(CMyVektor y, double x) {
        if(isHoehere) {
        CMyVektor ergebnisVek(y.getVektorDimension());
            for(int i = 0; i < y.getVektorDimension()-1; i++) {
                ergebnisVek.setVektorComponent(i,y.getVektorComponent(i+1));
            }
            ergebnisVek.setVektorComponent(y.getVektorDimension()-1,(*f_DGL_nterOrdnung)(y, x));
            return ergebnisVek;
        }
        else {
            return (*f_DGL_System)(y, x);
        }
    }

public:
    CDGLSolver(CMyVektor (*funktion)(CMyVektor y, double x))
        : f_DGL_System(funktion),f_DGL_nterOrdnung(nullptr), isHoehere(false) {}

    CDGLSolver(double (*funktion)(CMyVektor y, double x)) :
     f_DGL_System(nullptr), f_DGL_nterOrdnung(funktion), isHoehere(true){}

    // Getter für DGL-Typ
    bool getisHoehere() const {
        return isHoehere;
    }
    void setHoehere() {
        isHoehere = true;
    }
    void resetHoehere() {
        isHoehere = false;
    }
    int getSchritte() {return this->schritteHEun;}

    CMyVektor euler(double xStart, double xEnd, int nSteps, CMyVektor yStart) {
        double h = (xEnd - xStart) / nSteps;
        double x = xStart;
        CMyVektor y = yStart;

        // Iterationen
        for (int i = 0; i < nSteps; ++i) {
            CMyVektor dy = ableitungen(y, x);
            y = y + dy * h;
            x += h;
        }

        return y;
    }

    CMyVektor eulerAusgabeVariante(double xStart, double xEnd, int nSteps, CMyVektor yStart) {
        double h = (xEnd - xStart) / nSteps;
        double x = xStart;
        CMyVektor y = yStart;


        for (int i = 0; i < nSteps; ++i) {
            CMyVektor dy = ableitungen(y, x);
            y = y + dy * h;
            x += h;

            std::cout << "Schritt " << (i + 1) << ":" << std::endl;
            std::cout << "x = " << x << std::endl;

            std::cout << "y = (";
            for (int j = 0; j < y.getVektorDimension(); ++j) {
                std::cout << y.getVektorComponent(j);
                if (j + 1 < y.getVektorDimension()) std::cout << "; ";
            }
            std::cout << ")" << std::endl;

        }

        return y;
    }

    CMyVektor heun(double xStart, double xEnd, int nSteps, CMyVektor yStart) {
        double h = (xEnd - xStart) / nSteps;
        double x = xStart;
        CMyVektor y = yStart;

        // Schleife für alle Schritte
        for (int i = 0; i < nSteps; ++i) {
            CMyVektor dy_orig   = ableitungen(y, x);
            CMyVektor yTest     = y + dy_orig * h;
            CMyVektor dy_test   = ableitungen(yTest, x + h);
            CMyVektor dy_mittel = (dy_orig + dy_test) * 0.5;

            y = y + dy_mittel * h;
            x += h;
        }

        return y;
    }

    CMyVektor heunAusgabeVariante(double xStart, double xEnd, int nSteps, CMyVektor yStart) {
        double h = (xEnd - xStart) / nSteps;
        double x = xStart;
        CMyVektor y = yStart;

        for (int i = 0; i < nSteps; ++i) {
            // Heun-Schritte
            CMyVektor dy1 = ableitungen(y, x);
            CMyVektor yZwisch = y + dy1 * h;
            CMyVektor dy2 = ableitungen(yZwisch, x + h);
            CMyVektor dyMittel = (dy1 + dy2) * 0.5;

            y = y + dyMittel * h;
            x += h;

            // Einfache Ausgabe:
            std::cout << "Schritt " << (i + 1) << ":"<<std::endl;
            std::cout << "x = " << x << std::endl;
            std::cout << "y = (";

            // Diese Schleife gibt alle y-Komponenten aus
            for (int j = 0; j < y.getVektorDimension(); ++j) {
                std::cout << y.getVektorComponent(j);
                if (j + 1 < y.getVektorDimension()) std::cout << "; ";
            }

            std::cout << ")"<<std::endl;
        }

        return y;
    }





    void analysiereAbweichungDGL3(double xStart, double xEnd, CMyVektor yStart, double y_exakt) {
        std::vector<int> schritteListe = {10, 100, 1000, 10000};  // Schrittgrößen

        std::cout << "\nVergleich der Näherungen mit y(2) = " << y_exakt << "\n";
        std::cout << "---------------------------------------\n";

        // Schleife für jede Schrittgröße
        for (int i=0 ;i<schritteListe.size(); i++) {
            int n = schritteListe[i];

            CMyVektor eulerErgebnis = euler(xStart, xEnd, n, yStart);
            double abweichungEuler = eulerErgebnis.getVektorComponent(0) - y_exakt;
            std::cout << "Abweichung bei Euler mit " << n << " Schritten: " << abweichungEuler << "\n";


            CMyVektor heunErgebnis = heun(xStart, xEnd, n, yStart);
            double abweichungHeun = heunErgebnis.getVektorComponent(0) - y_exakt;
            std::cout << "Abweichung bei Heun mit " << n << " Schritten: " << abweichungHeun << "\n";

        }
    }



};