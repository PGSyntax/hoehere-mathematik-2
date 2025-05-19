//
// Created by legue on 06.04.2025.
//

#include <iostream>
#include<vector>
#include <cmath>


class CMyVektor {
private:

    unsigned int vektorDimension;
    std::vector<double> punkte;

public:

    CMyVektor(unsigned int dimension) : vektorDimension(dimension), punkte(vektorDimension,0.0) {
    }

    CMyVektor(const CMyVektor &v);

    CMyVektor& operator=(const CMyVektor &v);

    int getVektorDimension() const {

        return this->vektorDimension;
    }

    void setVektorComponent(int pos, double val) {
        if(pos>=vektorDimension || pos<0) {
            std::cout<<"Postition außerhalb der Dimension"<<std::endl;
            return;
        }
        punkte[pos] =val;
    }

    double getVektorComponent(int pos) const{
        if(pos>=vektorDimension || pos<0) {
            std::cout<<"Postition außerhalb der Dimension"<<std::endl;
            return 0.0;
        }
        return punkte[pos];
    }

    std::vector<double> getPunkte() const { return punkte; }

    double getVektorLength() const {
        double length = 0.0;
        for (int i = 0; i < punkte.size(); i++) {
            length+=(punkte[i]*punkte[i]);
        }
        return std::sqrt(length);
    }
    //-----------------------------------------------------------------------------------------------------------------------//

    CMyVektor gradient(CMyVektor x, double (*funktion)(CMyVektor x)) {
        double h = 1e-4; //10^-4
        CMyVektor gradient(x.getVektorDimension());
        for(int i = 0; i<x.getVektorDimension(); i++) {

            CMyVektor xiPlusH = x;
            xiPlusH.setVektorComponent(i,x.getVektorComponent(i)+h); //an stelle xi wird h addiert
            double diff = (funktion(xiPlusH) - funktion(x));
            gradient.setVektorComponent(i,diff/h);
        }
        return gradient;
    }

    //-----------------------------------------------------------------------------------------------------------------------//

CMyVektor gradientenVerfahren(CMyVektor startVektor, double (*funktion)(CMyVektor x),double lambda = 1.0) {
    CMyVektor x = startVektor;
    CMyVektor xNeu = x;

    // Abbruchkriterien
    double maxLaenge = 1E-5;
    int maxSchritte = 25;

    for (int counter = 1; counter <= maxSchritte; counter++) {
        std::cout<<"schritt "<<lambda <<" bei schritt :"<<counter<<std::endl;
        CMyVektor grad = gradient(x, funktion);
        if(grad.getVektorLength() < maxLaenge) break;
        xNeu = x + (lambda * grad);

        double fX = funktion(x);
        double fNeu = funktion(xNeu);

        if (fNeu > fX) {
            CMyVektor xNeu2 = x + (2 * lambda * grad);
            double fX2 = funktion(xNeu2);
            if(fX2>fNeu) {
                x=xNeu2;
                lambda *=2;
            }
            else {
                x = xNeu;
            }

        } else {
             do{
                lambda = lambda /2;
                xNeu = x + (lambda * grad);
                fNeu = funktion(xNeu);
            }while(fNeu < fX);
            x=xNeu;
        }
    }
    return x;
}


    friend CMyVektor operator+(const CMyVektor& a, const CMyVektor& b);
    friend CMyVektor operator*(double lambda, const CMyVektor& a);
    friend CMyVektor operator*(const CMyVektor &a, double lambda) {
        return lambda * a;
    }
    // Am Ende deiner Klassendefinition, innerhalb der Klasse:
    // Operator[] für den Zugriff auf die Elemente des Vektors

};


/*

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

CMyVektor::CMyVektor(const CMyVektor &v) : vektorDimension(v.getVektorDimension()), punkte(v.getPunkte()) {}*/






