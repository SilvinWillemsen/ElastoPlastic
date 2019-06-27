/*
 ==============================================================================
 
 ViolinString.h
 Created: 14 Nov 2018 3:10:58pm
 Author:  Silvin Willemsen
 
 ==============================================================================
 */

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"

enum InstrumentType
{
    twoStringViolin,
    sitar,
    bowedSitar,
    hurdyGurdy,
    dulcimer,
};

enum StringInterpolType
{
    noStringInterpol,
    linear,
    cubic,
};

enum BowModel
{
    exponential,
    elastoPlastic,
};


inline double exp1 (double x) {
    x = 1.0 + x / 1024;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x;
    return x;
}

using namespace std;
//==============================================================================
/*
 */
class ViolinString : public Component
{
public:
    ViolinString(double freq, double fs, int stringID, BowModel bowModel);
    ~ViolinString();
    
    void paint(Graphics &) override;
    void resized() override;
    
    void reset();
    void bow();
    void addJFc (double JFc, int index);
    
    void newtonRaphson();
    double getOutput (double ratio);
    
    void setFingerForce(double val) { ff = val; };
    void setFingerPosition(double fingerPos){ fp = clamp(fingerPos, 0, 1); };
    void setBow(bool val) { _isBowing.store(val); };
    void setVb(double val) { _Vb = val; };
    double getVb() { return _Vb; };
    double getQ() { return q; };
    void setFb(double val) { _Fb = val; };
    void setBowPos(double bpX, double bpY) { _bpX = bpX; _bpY = bpY; };
    void setFingerOn (bool val) { fingerOn = val; };
    
    double getOGFreq() { return ogFreq; };
    void setFrequency (double freq);
    
    double getNumPoints() { return N; };
    double getPrevStateAt (int idx) { return uPrev[idx]; };
    double getStateAt (int idx) { return u[idx]; };
    double getNextStateAt (int idx) { return uNext[idx]; };
    
    double getGridSpacing() { return h; };
    double getS0() { return s0; };
    
    void updateUVectors();
    
    void setRaisedCos (double exciterPos, double width, double force);
    bool isPicking() { return _isPicking;};
    
    Path generateStringPathAdvanced();
    
    void setConnection (int idx, double cp);
    int addConnection (double cp);
    
    int getCP (int idx) { return cpIdx[idx]; };
    int getCy (int idx) { return cy[idx]; };
    
    double clamp (double input, double min, double max);
    
    void mouseDown (const MouseEvent& e) override;
    void mouseDrag (const MouseEvent& e) override;
    void mouseUp (const MouseEvent& e) override;
    
    double linearInterpolation (double* uVec, int bowPosition, double alph);
    double cubicInterpolation (double* uVec, int bowPosition, double alph);
    
    int getStringID() { return stringID; };
    
    bool isStringBowing() { return _isBowing; };
    
    int sgn (double val) { return (0 < val) - (val < 0); };
    
    BowModel getModel() { return bowModel; };
    
    void setFn (double val) {
//        Fn = clamp(val, 0, 10);
        Fn = val;
        fC = mud * Fn;
        fS = mus * Fn;
    };
    double getFn() { return Fn; };
    double getFC() { return fC; };
    double getFS() { return fS; };
    double getSig0() { return sig0; };
    double getZ() { return z; };
    double getzBA() { return z_ba; };
    void scaleVisuals (double val) { visualScale = val; };
    
    void setNoise (double val) { sig3 = val; };
    
    void calcZDot();
    
private:
    double fs, freq, c, L, Eyoung, Iner, r, csA, rho, k, kHalf, s0, s1, B, kappa, h, N, muSq, kOh, gOh, a, BM, tol, q, qPrev, qPrevIt, b, eps, fp, B1, B2, b1, b2, A1, A1ss, A2, A3, A4, A5, D, E;
    double velCalcDiv = 0; //Velocity calculation division
    //Elasto plastic vars
    double z, zPrev, zPrevIt, zDot, zDotNext, zDotPrev, an, anPrev, scaleFact, fnl, z_ba, Fn, fC, fS, sig0, sig1, sig2, sig3, sig3w, oOSig0, E2, oOstrvSq, zss, zssNotAbs, oOZss, oOZssMinZba, dz_ss, dz_ssAbs, strv, espon, alpha, dalpha_v, dalpha_z, d_fnlv, d_fnlz, d_fnl, arg, mus, mud, K1, vRelTemp, zTemp, g1, g2, dg1v, dg1z, dg2v, dg2z, determ, invJac;
    double ff = 0.7;
    Random rand;
    atomic<double> lambdaSq;
    atomic<double> _Vb, _Fb, _bpX, _bpY;
    
    atomic<double> bowPos;
    double ogFreq; //original frequency
    vector<int> cx;
    vector<int> cy;
    vector<int> cpIdx;
    atomic<bool> _isBowing{false};
    atomic<bool> _isPicking{false};
    
    double* uNext;
    double* u;
    double* uPrev;
    int uNextPtrIdx = 0;
    vector<vector<double>> uVecs;
    
    //NR interpolation variables
    double uI = 0;
    double uIPrev = 0;
    double uI1 = 0;
    double uI2 = 0;
    double uIM1 = 0;
    double uIM2 = 0;
    double uIPrev1 = 0;
    double uIPrevM1 = 0;
    
    bool active = false;
    
    unsigned long count;
    unsigned long t = 0;
    
    bool fingerOn = false;
    int fpx = 0;
    
    StringInterpolType interpolation = noStringInterpol;
    int cpMoveIdx = -1;
    int cpMR = 10; //connection point move range
    
    int stringID;
    
    BowModel bowModelInit;
    BowModel bowModel;
    
    double visualScale = 5.0; 
    
    double oneOverSix = 1/6.0;
    long limitCount = 0;
    int maxI = 0;
    
    bool simplySupported = true;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ViolinString)
};


