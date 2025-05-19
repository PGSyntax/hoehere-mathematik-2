#include <iostream>
#include <cmath>
#include "C_DGLSolver.cpp"



CMyVektor dglSystem(CMyVektor y, double x) {
    CMyVektor dy(2);
    dy.setVektorComponent(0, 2 * y.getVektorComponent(1) - x * y.getVektorComponent(0));
    dy.setVektorComponent(1, y.getVektorComponent(0) * y.getVektorComponent(1) - 2*x*x*x);
    return dy;
}

double dglDritteOrdnung(CMyVektor y, double x) {
    double y_ = y.getVektorComponent(0);   // y
    double y1 = y.getVektorComponent(1);   // y'
    double y2 = y.getVektorComponent(2);   // y''
    return 2 * x * y1 * y2 + 2 * y_ * y_ * y1;
}


int main() {
    CDGLSolver solverDGL1(dglSystem);  // System erster Ordnung
    CDGLSolver solverDGL3(dglDritteOrdnung);  // DGL dritter Ordnung

    CMyVektor y0(2);
    y0.setVektorComponent(0, 0.0);
    y0.setVektorComponent(1, 1.0);

    std::cout << "Wählen Sie das Verfahren (1 = Euler, 2 = Heun, 3 = DGL 3. Ordnung mit Abweichung): ";
    int choice;
    std::cin >> choice;
    std::cout << "\n";

    if (choice == 1) {
        std::cout << "=== Euler-Verfahren ===\n";
        solverDGL1.eulerAusgabeVariante(0.0, 2.0, 100, y0);
    } else if (choice == 2) {
        std::cout << "=== Heun-Verfahren ===\n";
        solverDGL1.heunAusgabeVariante(0.0, 2.0, 100, y0);
    } else if (choice == 3) {
        std::cout << "=== DGL dritter Ordnung: Abweichungsanalyse ===\n";

        // Startwerte für die DGL dritter Ordnung
        CMyVektor yStart(3);
        yStart.setVektorComponent(0, 1.0);   // y(1)
        yStart.setVektorComponent(1, -1.0);  // y'(1)
        yStart.setVektorComponent(2, 2.0);   // y''(1)

        // Exakte Lösung für y(2) ist 0.5
        solverDGL3.analysiereAbweichungDGL3(1.0, 2.0, yStart, 0.5);
    }


    return 0;
}


/*
int main() {
    // DGL-System definieren (zwei Gleichungen)
    CMyVektor (*fSystem)(CMyVektor, double) = [](CMyVektor y, double x) {
        CMyVektor dy(2);
        dy.setVektorComponent(0, 2 * y.getVektorComponent(1) - x * y.getVektorComponent(0));
        dy.setVektorComponent(1, y.getVektorComponent(0) * y.getVektorComponent(1) - 2*x*x*x);
        return dy;
    };

    CDGLSolver solver(fSystem);
    CMyVektor y0(2);
    y0.setVektorComponent(0, 0.0);
    y0.setVektorComponent(1, 1.0);
    int nSteps = 100;

    std::cout << "Wählen Sie das Verfahren (1 = Euler, 2 = Heun): ";
    int choice;
    std::cin >> choice;
    std::cout << "\n";

    if (choice == 1) {
        std::cout << "=== Euler-Verfahren ===\n";
        solver.euler(0.0, 2.0, nSteps, y0);
    } else if (choice == 2) {
        std::cout << "=== Heun-Verfahren ===\n";
        solver.heun(0.0, 2.0, nSteps, y0);
    }

    return 0;
}


*/





/*

// ... hier deine Klasse CDGLSolver wie gehabt ...
int main() {
    // Beispiel: DGL-System y'1=2y2 - x y1, y'2 = y1*y2 -2x^3
    auto fSystem = [](CMyVektor y, double x) {
        CMyVektor dy(2);
        dy.setVektorComponent(0, 2 * y.getVektorComponent(1) - x * y.getVektorComponent(0));
        dy.setVektorComponent(1, y.getVektorComponent(0) * y.getVektorComponent(1) - 2*x*x*x);
        return dy;
    };
    CDGLSolver solverSys(fSystem);
    CMyVektor y0_sys(2);
    y0_sys.setVektorComponent(0, 0.0);
    y0_sys.setVektorComponent(1, 1.0);
    int nSteps = 100;

    std::cout << "=== Euler-Verfahren ===\n";
    solverSys.euler(0.0, 2.0, nSteps, y0_sys);

    std::cout << "=== Heun-Verfahren ===\n";
    solverSys.heun (0.0, 2.0, nSteps, y0_sys);

    return 0;
}

*/









/* p2
// Testfunktion für Aufgabe 2: f: R^4 -> R^3
CMyVektor testfunktion1(CMyVektor x) {
    CMyVektor ergebnis(3);
    ergebnis.setVektorComponent(0, x.getVektorComponent(0) * x.getVektorComponent(1) * exp(x.getVektorComponent(2)));
    ergebnis.setVektorComponent(1, x.getVektorComponent(1) * x.getVektorComponent(2) * x.getVektorComponent(3));
    ergebnis.setVektorComponent(2, x.getVektorComponent(3));
    return ergebnis;
}

// Testfunktion für Newtonverfahren: f: R^2 -> R^2
CMyVektor testfunktion2(CMyVektor x) {
    CMyVektor ergebnis(2);
    ergebnis.setVektorComponent(0, pow(x.getVektorComponent(0), 3) * pow(x.getVektorComponent(1), 3) - 2 * x.getVektorComponent(1));
    ergebnis.setVektorComponent(1, x.getVektorComponent(0) - 2);
    return ergebnis;
}

int main() {

    std::cout << "\nBerechnung der Jacobi-Matrix fuer Funktion f: R^4 -> R^3:" << std::endl;

    CMyVektor x1(4);
    x1.setVektorComponent(0, 1);
    x1.setVektorComponent(1, 2);
    x1.setVektorComponent(2, 0);
    x1.setVektorComponent(3, 3);

    //3x4 Matrix
    CMyMatrix MtxA2(3, 4);
   // MtxA2.print();

    CMyMatrix jacobiMatrix = MtxA2.jacobi(x1, testfunktion1);

    // Ausgabe der Jacobi-Matrix
    std::cout << "Jacobi-Matrix Jf(x) =\n";
    for (int i = 0; i < jacobiMatrix.getDimension(); ++i) {
        for (int j = 0; j < jacobiMatrix.getSpalten(); ++j) {
            std::cout << jacobiMatrix.getMatrixComponent(i, j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nStarte Newtonverfahren fuer Funktion f: R^2 -> R^2:" << std::endl;

    // Startwert für Newton-Verfahren
    CMyVektor startwert(2);
    startwert.setVektorComponent(0, 1);
    startwert.setVektorComponent(1, 1);
    //newton mtx
    CMyMatrix MtxA3(2, 2);
    //MtxA3.print(); //zeilenausgabe für debug

    CMyVektor ergebnis = MtxA3.newtonVerfahren(startwert, testfunktion2);

    std::cout << "\nGefundene Nullstelle: ("
              << ergebnis.getVektorComponent(0) << ", "
              << ergebnis.getVektorComponent(1) << ")\n";

    return 0;
}

/*















/*
// Funktion f(x, y) = sin(xy) + sin(x) + cos(y)
double f(CMyVektor x) {

    double x1 = x.getVektorComponent(0);
    double x2 = x.getVektorComponent(1);
    return sin(x1 * x2) + sin(x1) + cos(x2);
}

// Funktion g(x1, x2, x3) = −(2x1^2 − 2x1x2 + x2^2 + x3^2 − 2x1 − 4x3)
double g(CMyVektor x) {

    double x1 = x.getVektorComponent(0);
    double x2 = x.getVektorComponent(1);
    double x3 = x.getVektorComponent(2);
    return -(2 * x1 * x1 - 2 * x1 * x2 + x2 * x2 + x3 * x3 - 2 * x1 - 4 * x3);
}

int main() {
    std::cout << "Test 1: Funktion f(x, y) = sin(xy) + sin(x) + cos(y)" << std::endl;
    CMyVektor x1(2);
    x1.setVektorComponent(0, 0.2);
    x1.setVektorComponent(1, -2.1);

    CMyVektor result1 = x1.gradientenVerfahren(x1, f);
    std::cout << "Ergebnis nach Gradientenverfahren für f(x, y):" << std::endl;
    for (int i = 0; i < result1.getVektorDimension(); ++i) {
        std::cout << "x[" << i << "] = " << result1.getVektorComponent(i) << std::endl;
    }

    std::cout << "\nTest 2: Funktion g(x1, x2, x3) = −(2x1^2 − 2x1x2 + x2^2 + x3^2 − 2x1 − 4x3)" << std::endl;
    CMyVektor x2(3);
    x2.setVektorComponent(0, 0.0);
    x2.setVektorComponent(1, 0.0);
    x2.setVektorComponent(2, 0.0);

    CMyVektor result2 = x2.gradientenVerfahren(x2, g, 0.1);
    std::cout << "Ergebnis nach Gradientenverfahren für g(x1, x2, x3):" << std::endl;
    for (int i = 0; i < result2.getVektorDimension(); ++i) {
        std::cout << "x[" << i << "] = " << result2.getVektorComponent(i) << std::endl;
    }

    return 0;
}
*/

CMyVektor operator+(const CMyVektor& a, const CMyVektor& b) {
    CMyVektor neuerVektor(a.getVektorDimension());

    for(int i = 0; i < a.getVektorDimension(); i++) {
        neuerVektor.punkte[i]= b.punkte[i]+a.punkte[i];
    }
    return neuerVektor;
}

CMyVektor operator*(double lambda, const CMyVektor& a) {
    if(a.getVektorDimension()>0) {
        CMyVektor neuerVektor(a.getVektorDimension());

        for(int i=0;i<a.getVektorDimension();i++) {
            neuerVektor.punkte[i] = (lambda*a.punkte[i]);
        }
        return neuerVektor;
    }
    else {
        std::cout<<"eingabe vektor ist leer"<<std::endl;
        return CMyVektor(0);
    }

}

CMyVektor& CMyVektor::operator=(const CMyVektor& a) {
    this->punkte = a.getPunkte();
    return *this;
}

CMyVektor::CMyVektor(const CMyVektor &v) : vektorDimension(v.getVektorDimension()), punkte(v.getPunkte()) {}

