/*
 ==============================================================================
 
 ViolinString.cpp
 Created: 14 Nov 2018 3:10:58pm
 Author:  Silvin Willemsen
 
 ==============================================================================
 */

#include "../JuceLibraryCode/JuceHeader.h"
#include "ViolinString.h"

//==============================================================================
ViolinString::ViolinString (double freq, double fs, int stringID, BowModel bowModelInit) : fs(fs), ogFreq(freq), stringID(stringID), bowModelInit (bowModelInit)
{
    uVecs.resize(3);
    bowModel = bowModelInit;
    _bpX = 0.5;
    _bpY = 0;
    
    c = ogFreq * 2; // Wave speed
    k = 1 / fs;       // Time-step
    kHalf = 0.5 * k;
    L = 1;
    
    rho = 7850;
    r = 0.0005;
    csA = double_Pi * r * r;
    Iner = double_Pi * r * r * r * r / 4.0;
    Eyoung = 2e11;
    
    double T = c * c * rho * csA;
    
    s0 = 1 * rho * csA;     // Frequency-independent damping
    s1 = 0.005 * rho * csA; // Frequency-dependent damping

    kappa = sqrt(Eyoung * Iner / (rho * csA));
    // Grid spacing
    //    if (stringType == bowedString && instrumentType == twoStringViolin)
    //    {
    //        N = 10;
    //    }
    //    else
    
    h = sqrt((c * c * k * k + 4.0 * s1 * k + sqrt(pow(c * c * k * k + 4.0 * s1 * k, 2.0) + 16.0 * kappa * kappa * k * k)) * 0.5);
    N = floor(L / h); // Number of gridpoints
    
    h = L / N; // Recalculate gridspacing
    
    // Initialise vectors
    //    vector<double> dummyVector (N, 0);
    for (int i = 0; i < 3; ++i)
    {
        uVecs[i].resize(N);
    }
    
    uNext = &uVecs[uNextPtrIdx][0];
    u = &uVecs[1][0];
    uPrev = &uVecs[2][0];
    
    // Courant numbers
    lambdaSq = c * c / (h * h);
    muSq = kappa * kappa / (h * h * h * h);
    
    kOh = (kappa * kappa) / (h * h * h * h);
    gOh = (c * c) / (h * h);
    
    // Simple (exponential) Bow Model
    a = 100; // Free parameter
    BM = sqrt(2 * a) * exp1(0.5);
    
    _Vb = 0.1; // Bowing speed
    _Fb = 80;  // Bowing force / total mass of bow;
    
    // Elasto-Plastic bow model
    
    //// the Contact Force (be with you) //////
    mus = 0.8; // static friction coeff
    mud = 0.3; // dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
    strv = 0.1;      // "stribeck" velocity
    
    Fn = 1;    // Normal force
    
    fC = mud * Fn; // coulomb force
    fS = mus * Fn; // stiction force
    
    sig0 = 10000;                   // bristle stiffness
    sig1 = 0.001*sqrt(sig0);          // bristle damping
    sig2 = 0.4;                     // viscous friction term
    sig3 = 0;                       // noise term
    oOstrvSq = 1 / (strv * strv);   // One over strv^2
    z_ba = 0.7 * fC / sig0;         // break-away displacement (has to be < f_c/sigma_0!!)
    
    // Initialise variables for Newton Raphson
    tol = 1e-7;
    qPrev = 0;
    zPrev = 0;
    zDotPrev = 0;
    anPrev = 0;
    fp = 0;
    
    // FDS precalculations
    B1 = s0 / k;
    B2 = (2 * s1) / (k * h * h);
    
    b1 = 2.0 / (k * k);
    b2 = (2 * s1) / (k * h * h);
    
    D = 1.0 / (rho * csA / (k * k) + s0 / k);
    
    A1 = 2 * rho * csA / (k * k) - 2 * T / (h * h) - 6 * Eyoung * Iner / (h * h * h * h) - 2 * B2;
    A1ss = 2 * rho * csA / (k * k) - 2 * T / (h * h) - 5 * Eyoung * Iner / (h * h * h * h) - 2 * B2;
//    A1 = (2 - 2 * lambdaSq - 6 * muSq - 2 * B2) * (rho * csA / (k * k)) / (rho * csA / (k * k) + s0 / k);
    A2 = T / (h * h) + 4 * Eyoung * Iner / (h * h * h * h) + B2;
    A3 = Eyoung * Iner / (h * h * h * h);
    A4 = B1 - (rho * csA) / (k * k) + 2 * B2;
    A5 = B2;
    
    A1 *= D;
    A1ss *= D;
    A2 *= D;
    A3 *= D;
    A4 *= D;
    A5 *= D;
    
    E = k * k * (1 / h) * BM;
    E2 = 1.0 / (h * (rho * csA / (k * k) + s0 / k));
    
    reset();
    q = 0;
    
    velCalcDiv =  1 / (sig2 + scaleFact * (2/k + 2 * s0));
    oOSig0 = 1 / sig0;
    _isBowing =  true;
    std::cout << N << std::endl;
}

void ViolinString::reset()
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < N; ++j)
            uVecs[i][j] = 0.0;
    qPrev = 0;
    q = 0;
    zPrev = 0;
    zDotPrev = 0;
    anPrev = 0;
    zPrevIt = 0;
    qPrevIt = 0;
    z = 0;
    zDot = 0;
    std::cout << "RESET" << std::endl;
}

ViolinString::~ViolinString()
{
}

void ViolinString::paint(Graphics &g)
{
    /* This demo code just fills the component's background and
     draws some placeholder text to get you started.
     
     You should replace everything in this method with your own
     drawing code..
     */
    bowModel = bowModelInit;
    g.setColour(bowModel == elastoPlastic ? Colours::cyan : Colours::limegreen);
    g.strokePath(generateStringPathAdvanced(), PathStrokeType(2.0f));

    g.setColour(Colours::orange);
    for (int c = 0; c < cpIdx.size(); ++c)
    {
        g.drawEllipse(floor(cpIdx[c] * getWidth() / N - 5), floor(cy[c] - 5), 10, 10, 2);
    }
    
    g.setColour(Colours::yellow);
    g.fillEllipse(fp * getWidth() - 5, getHeight() / 2.0 - 5, 10, 10);
    
    // draw bow
    Colour c = Colours::yellow;
    float alph = static_cast<float>((Fn + 20) / 80.0);
//    std::cout << alph << std::endl;
    g.setColour (c.withAlpha (alph));
//    double opa = 90.0 / 100.0;
//    if (opa >= 1.0)
//    {
//        g.setOpacity(1.0);
//    }
//    else
//    {
//        g.setOpacity(opa);
//    }
    g.fillRect(floor(_bpX.load() * getWidth()), floor(_bpY.load() * getHeight()) - getHeight() / 2.0, 10, getHeight());
    g.setColour(Colour::greyLevel(0.5f).withAlpha(0.5f));
    for (double i = -12.0; i < 12.0; ++i)
    {
        double val = (1 - (pow(2.0, (i / 12.0)) - 1)) * getWidth() * 0.5;
        //            std::cout << val << std::endl;
        
        g.drawLine(val, 0, val, getHeight(), 2);
    }
    g.setColour(Colour::fromRGBA(255, 255, 0, 127));
    g.drawLine(getWidth() / 2.0, 0, getWidth() / 2.0, getHeight(), 2);
    
}

void ViolinString::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..
}

void ViolinString::bow()
{
    std::vector<double> uSave (N, 0);
    for (int l = 0; l < N; ++l)
    {
        uSave[l] = u[l];
    }
    
    double Fb = _Fb.load();
    bowPos.store(clamp(_bpX.load() * N, 2, N - 3));
    int bp = floor(bowPos.load());
    bool isBowing = _isBowing;
    
    if (isBowing)
    {
//        sig3w = (rand.nextFloat() * 2 - 1) * sig3;
        sig3w = 0;
        newtonRaphson();
    }
    //    if (pluck && pluckIdx < pluckLength)
    //    {
    //
    //    }
    double excitation;
    if (bowModel == exponential)
    {
        excitation = E * Fb * q * exp1(-a * q * q);
    }
    else if (bowModel == elastoPlastic && isBowing)
    {
        excitation = E2 * (sig0 * z + sig1 * zDot + sig2 * q + sig3w); //* (rho * csA);
    } else {
        excitation = 0;
    }
    
    for (int l = 2; l < N - 2; ++l)
    {
        uNext[l] = A1 * u[l] + A2 * (u[l + 1] + u[l - 1]) - A3 * (u[l + 2] + u[l - 2]) + A4 * uPrev[l] - A5 * (uPrev[l + 1] + uPrev[l - 1]);
    }
    
    if (simplySupported)
    {
        int l = 1;
        uNext[l] = A1ss * u[l] + A2 * (u[l + 1] + u[l - 1]) - A3 * (u[l + 2]) + A4 * uPrev[l] - A5 * (uPrev[l + 1] + uPrev[l - 1]);
        l = N - 2;
        uNext[l] = A1ss * u[l] + A2 * (u[l + 1] + u[l - 1]) - A3 * (u[l - 2]) + A4 * uPrev[l] - A5 * (uPrev[l + 1] + uPrev[l - 1]);
    }
    if (isBowing)
    {
        double alpha = bowPos - floor(bowPos);

        if (interpolation == noStringInterpol)
        {
            uNext[bp] = uNext[bp] - excitation;
        }
        else if (interpolation == linear)
        {
            uNext[bp] = uNext[bp] - excitation * (1-alpha);
            
            if (bp < N - (simplySupported ? 2 : 3))
                uNext[bp + 1] = uNext[bp + 1] - excitation * alpha;
        }
        else if (interpolation == cubic)
        {
            if (bp > (simplySupported ? 2 : 3))
                uNext[bp - 1] = uNext[bp - 1] - excitation * (alpha * (alpha - 1) * (alpha - 2)) * -oneOverSix;
            
            uNext[bp] = uNext[bp] - excitation * ((alpha - 1) * (alpha + 1) * (alpha - 2)) * 0.5;
            
            if (bp < N - (simplySupported ? 2 : 3))
                uNext[bp + 1] = uNext[bp + 1] - excitation * (alpha * (alpha + 1) * (alpha - 2)) * -0.5;
            if (bp < N - (simplySupported ? 3 : 4))
                uNext[bp + 2] = uNext[bp + 2] - excitation * (alpha * (alpha + 1) * (alpha - 1)) * oneOverSix;
        }
    }
   
    if (ff > 1)
    {
        ff = 1;
        std::cout << "wait" << std::endl;
    }
    ///////// heuristic interpolation
    int fingerPos = floor(fp * N - 1);
    double scale = 1;
    double alphaFP = (fp * N - 1 - fingerPos) * scale;
    if (fingerPos > 1)
    {
        uNext[fingerPos - 1] *= 0;
    }
    if (fingerPos < N - 2)
    {
        uNext[fingerPos] *= 0;
    }
    if (fingerPos < N - 3)
    {
        uNext[fingerPos + 1] = uNext[fingerPos + 1] * (1 - alphaFP);
    }
    if (isnan(uNext[static_cast<int>(N*0.5)]))
    {
        reset();
        return;
    }
}

void ViolinString::newtonRaphson()
{
//    double Vb = _Vb.load();
    double Vb = _Vb.load();
    double Fb = _Fb.load();
    int bp = floor(bowPos.load());
    double alpha = bowPos.load() - bp;
    
    if (interpolation == noStringInterpol)
    {
        uI = u[bp];
        uIPrev = uPrev[bp];
        uI1 = u[bp + 1];
        uI2 = u[bp + 2];
        uIM1 = u[bp - 1];
        uIM2 = u[bp - 2];
        uIPrev1 = uPrev[bp + 1];
        uIPrevM1 = uPrev[bp - 1];
    }
    else if (interpolation == linear)
    {
        uI = linearInterpolation(u, bp, alpha);
        uIPrev = linearInterpolation(uPrev, bp, alpha);
        uI1 = linearInterpolation(u, bp + 1, alpha);
        uI2 = linearInterpolation(u, bp + 2, alpha);
        uIM1 = linearInterpolation(u, bp - 1, alpha);
        uIM2 = linearInterpolation(u, bp - 2, alpha);
        uIPrev1 = linearInterpolation(uPrev, bp + 1, alpha);
        uIPrevM1 = linearInterpolation(uPrev, bp - 1, alpha);
    }
    else if (interpolation == cubic)
    {
        uI = cubicInterpolation(u, bp, alpha);
        uIPrev = cubicInterpolation(uPrev, bp, alpha);
        uI1 = cubicInterpolation(u, bp + 1, alpha);
        uI2 = cubicInterpolation(u, bp + 2, alpha);
        uIM1 = cubicInterpolation(u, bp - 1, alpha);
        uIM2 = cubicInterpolation(u, bp - 2, alpha);
        uIPrev1 = cubicInterpolation(uPrev, bp + 1, alpha);
        uIPrevM1 = cubicInterpolation(uPrev, bp - 1, alpha);
    }
    
//    b = 2.0 / k * Vb - b1 * (uI - uIPrev) - gOh * (uI1 - 2 * uI + uIM1) + kOh * (uI2 - 4 * uI1 + 6 * uI - 4 * uIM1 + uIM2) + 2 * s0 * Vb - b2 * ((uI1 - 2 * uI + uIM1) - (uIPrev1 - 2 * uIPrev + uIPrevM1)) - sig3w;
    eps = 1;
    int i = 0;
//    std::cout << fS << std::endl;
    if (bowModel == exponential)
    {
        b = 2.0 / k * Vb - b1 * (uI - uIPrev) - gOh * (uI1 - 2 * uI + uIM1) + kOh * (uI2 - 4 * uI1 + 6 * uI - 4 * uIM1 + uIM2) + 2 * s0 * Vb - b2 * ((uI1 - 2 * uI + uIM1) - (uIPrev1 - 2 * uIPrev + uIPrevM1));
        
//        std::cout<< "Exponential model" << std::endl;
        while (eps > tol)
        {
            q = qPrev - (Fb * BM * qPrev * exp1(-a * qPrev * qPrev) + 2 * qPrev / k + 2 * s0 * qPrev + b) / (Fb * BM * (1 - 2 * a * qPrev * qPrev) * exp1(-a * qPrev * qPrev) + 2 / k + 2 * s0);
            eps = std::abs(q - qPrev);
            qPrev = q;
            ++i;
            if (i > 10000)
            {
                std::cout << "Nope" << std::endl;
            }
        }
    }
    else if (bowModel == elastoPlastic)
    {
        b = 2.0 / k * Vb - b1 * (uI - uIPrev) - gOh * (uI1 - 2 * uI + uIM1) + kOh * (uI2 - 4 * uI1 + 6 * uI - 4 * uIM1 + uIM2) + 2 * s0 * Vb - b2 * ((uI1 - 2 * uI + uIM1) - (uIPrev1 - 2 * uIPrev + uIPrevM1));
        z_ba = 0.7 * fC * oOSig0;
        // b
        while (eps > tol && i < 50 && fC > 0)
        {
            calcZDot();
            
            g1 = (2.0 / k + 2 * s0) * q + (sig0 * z + sig1 * zDot + sig2 * q + sig3w) / (rho * csA * h) + b;
            g2 = zDot - an;
            
            // compute derivatives
            
            // dz_ss/dv
            dz_ss = (-2 * abs(q) * oOstrvSq * oOSig0) * (fS-fC) * espon;
            dz_ssAbs = sgn(zss) * dz_ss;
            
            dalpha_v = 0; //d(alpha)/dv
            dalpha_z = 0; //d(alpha)/dz
            zss = abs(zss);
            if ((sgn(z)==sgn(q)) && (abs(z)>z_ba) && (abs(z)<zss) )
            {
                double cosarg = cos(sgn(z) * arg);
                dalpha_v = 0.5 * double_Pi * cosarg * dz_ssAbs * (z_ba - abs(z)) * oOZssMinZba * oOZssMinZba;
                dalpha_z = 0.5 * double_Pi * cosarg * sgn(z) * oOZssMinZba;
            }
            zss = zssNotAbs;
            d_fnlv = 1 - z * ((alpha + q * dalpha_v) * zss - dz_ss * alpha * q) * oOZss * oOZss;
            d_fnlz = -q * oOZss * (z * dalpha_z + alpha);
//            d_fnl = d_fnlv * K1 + d_fnlz * kHalf;
            
            dg1v = 2.0 / k + 2 * s0 + sig1 / (rho * csA * h) * d_fnlv + sig2 / (rho * csA * h);
            dg1z = sig0 / (rho * csA * h) + sig1 / (rho * csA * h) * d_fnlz;
            dg2v = d_fnlv;
            dg2z = d_fnlz - 2.0 / k;
            
            determ = dg1v * dg2z - dg1z * dg2v;
            qPrevIt = q;
            zPrevIt = z;
            q = q - (1 / determ) * (dg2z * g1 - dg1z * g2);
            z = z - (1 / determ) * (-dg2v * g1 + dg1v * g2);
            
            eps = sqrt((q-qPrevIt)*(q-qPrevIt) + (z-zPrevIt)*(z-zPrevIt));
            i = i + 1;
        }
        if (i == 50)
        {
            ++limitCount;
            std::cout << Fn << " Limit! " << limitCount <<  std::endl;
        }
        
        calcZDot();
        
        zPrev = z;
        zDotPrev = zDot;
        anPrev = an;
    }
}

void ViolinString::calcZDot()
{
    espon = exp1 (-((q * q) * oOstrvSq));         //exponential function
    zss = sgn(q) * (fC + (fS - fC) * espon) * oOSig0;   //steady state curve: z_ss(v)
    //            std::cout << zss << std::endl;
    if (q==0)
        zss = fS * oOSig0;
    
    // elasto-plastic function \alpha (v,z)
    alpha=0;
    
    oOZss = 1 / zss; // should use the absolute zss
    zssNotAbs = zss;
    zss = abs(zss);
    
    oOZssMinZba = 1 / (zss - z_ba); // should use the absolute zss
    
    if (sgn(z)==sgn(q))
    {
        if ((abs(z)>z_ba) && (abs(z)<zss))
        {
            arg = double_Pi * (z - sgn(z) * 0.5 * (zss + z_ba)) * oOZssMinZba;
            alpha = 0.5 * (1 + sin(sgn(z) * arg));
        }
        else if (abs(z)>=zss)
        {
            alpha=1;
        }
    }
    zss = zssNotAbs;
    an = 2.0 / k * (z - zPrev) - anPrev;
    
    // non-linear function estimate
    zDot = q * (1 - alpha * z * oOZss);
}

void ViolinString::addJFc(double JFc, int index)
{
    uNext[index] = uNext[index] + JFc;
}

double ViolinString::getOutput(double ratio)
{
    int index = floor(ratio * N);
    return uNext[index];
}

void ViolinString::updateUVectors()
{
    uPrev = u;
    u = uNext;
    
    uNextPtrIdx = (uNextPtrIdx + 2) % 3;
    uNext = &uVecs[uNextPtrIdx][0];
}

Path ViolinString::generateStringPathAdvanced()
{
    auto stringBounds = getHeight() / 2.0;
    Path stringPath;
    stringPath.startNewSubPath(0, stringBounds);
    
    auto spacing = getWidth() / static_cast<double>(N);
    auto x = spacing;
    
    for (int y = 1; y < N-1; y++)
    {
        int visualScaling = (bowModel == elastoPlastic ? 100000 : 10000) * visualScale;
        float newY = uNext[y] * visualScaling + stringBounds;
        if (isnan(newY))
            newY = 0;
        stringPath.lineTo(x, newY);
        for (int c = 0; c < cpIdx.size(); ++c)
        {
            if (y == cpIdx[c])
            {
                cx[c] = x;
                cy[c] = newY;
            }
        }
        if (y == floor(fp * N))
        {
            fpx = x;
        }
        x += spacing;
    }
    stringPath.lineTo(getWidth(), stringBounds);
    
    return stringPath;
}

void ViolinString::setFrequency (double freq)
{
    c = 2 * freq;
    lambdaSq = (c * c * k * k) / (h * h);
    A1 = 2 - 2 * lambdaSq - 6 * muSq - 2 * B2;
    A2 = lambdaSq + 4 * muSq + B2;
}

double ViolinString::clamp (double input, double min, double max)
{
    if (input > max)
        return max;
    else if (input < min)
        return min;
    else
        return input;
}

void ViolinString::mouseDown(const MouseEvent &e)
{
//    std::cout << (bowModel == exponential ? "Exponential" : "Elasto-Plastic") << std::endl;
    if (ModifierKeys::getCurrentModifiers() == ModifierKeys::leftButtonModifier)
    {
//        _Vb = 0.1;
//        _Fb = 80;
        setBow(true);
    }
    if (e.y >= (getHeight() / 2.0) - cpMR && e.y <= (getHeight() / 2.0) + cpMR && ModifierKeys::getCurrentModifiers() == ModifierKeys::altModifier + ModifierKeys::leftButtonModifier)
    {
        bool cpMove = false;
        for (int i = 0; i < cpIdx.size(); ++i)
        {
            if (getWidth() * cpIdx[i] / static_cast<double>(N) >= e.x - cpMR && getWidth() * cpIdx[i] / static_cast<double>(N) <= e.x + cpMR)
            {
                cpMove = true;
                cpMoveIdx = i;
                break;
            }
        }
        // if the location of the click didn't contain an existing connection, create a new one
        if (!cpMove)
        {
            //create new connection
            std::cout << "check" << std::endl;
        }
    }
}

void ViolinString::mouseDrag(const MouseEvent &e)
{
    double maxVb = 0.2;
    
    if (cpMoveIdx != -1 || ModifierKeys::getCurrentModifiers() == ModifierKeys::altModifier + ModifierKeys::leftButtonModifier)
    {
        if (cpMoveIdx != -1)
        {
            double cp = e.x <= 0 ? 0 : (e.x < getWidth() ? e.x / static_cast<double>(getWidth()) : 1);
            cpIdx[cpMoveIdx] = floor(cp * N);
        }
    }
    else if (ModifierKeys::getCurrentModifiers() == ModifierKeys::ctrlModifier + ModifierKeys::leftButtonModifier)
    {
        fp = e.x <= 0 ? 0 : (e.x < getWidth() ? e.x / static_cast<double>(getWidth()) : 1);
    }
    
    else if (cpMoveIdx == -1)
    {
        float bowVelocity = (e.y - getHeight() * 0.5) / (static_cast<double>(getHeight())) * maxVb * 2;
        setVb (bowVelocity);
        
        float bowPositionX = e.x <= 0 ? 0 : (e.x < getWidth() ? e.x / static_cast<double>(getWidth()) : 1);
        float bowPositionY = e.y <= 0 ? 0 : (e.y >= getHeight() ? 1 : e.y / static_cast<double>(getHeight()));
        
        _bpX.store(bowPositionX);
        _bpY.store(bowPositionY);
    }
}

void ViolinString::mouseUp(const MouseEvent &e)
{
    setBow(false);
    cpMoveIdx = -1;
}

double ViolinString::linearInterpolation(double* uVec, int bp, double alph)
{
    return uVec[bp] * (1 - alph) + uVec[bp + 1] * alph;
}

double ViolinString::cubicInterpolation(double* uVec, int bp, double alph)
{
//    double val1 = uVec[bp - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0
//    + uVec[bp] * ((alph - 1) * (alph + 1) * (alph - 2)) / 2.0
//    + uVec[bp + 1] * (alph * (alph + 1) * (alph - 2)) / -2.0
//    + uVec[bp + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
//
    double val = 0;
    if (bp > simplySupported ? 2 : 3)
        val = val + uVec[bp - 1] * (alph * (alph - 1) * (alph - 2)) * -oneOverSix;
    
    val = val + uVec[bp] * ((alph - 1) * (alph + 1) * (alph - 2)) * 0.5;
    
    if (bp < N - (simplySupported ? 2 : 3))
        val = val + uVec[bp + 1] * (alph * (alph + 1) * (alph - 2)) * -0.5;
    if (bp < N - (simplySupported ? 3 : 4))
        val = val + uVec[bp + 2] * (alph * (alph + 1) * (alph - 1)) * oneOverSix;
//    std::cout << val - val1 << std::endl;
    return val;
}
