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
ViolinString::ViolinString (double freq, double fs, int stringID, BowModel bowModel) : fs(fs), ogFreq(freq), stringID(stringID), bowModelInit(bowModel)
{
    // In your constructor, you should add any child components, and
    // initialise any special settings that your component needs.
    uVecs.resize(3);
    
    _bpX = 0.5;
    _bpY = 0;
    
    c = ogFreq * 2; // Wave speed
    k = 1 / fs;       // Time-step
    
    s0 = 1;     // Frequency-independent damping
    s1 = 0.005; // Frequency-dependent damping
    
    //    B = 0.0001;                             // Inharmonicity coefficient
    //    kappa = sqrt (B) * (gamma / double_Pi); // Stiffness Factor
    kappa = 2.0;
    // Grid spacing
    //    if (stringType == bowedString && instrumentType == twoStringViolin)
    //    {
    //        N = 10;
    //    }
    //    else
    
    h = sqrt((c * c * k * k + 4.0 * s1 * k + sqrt(pow(c * c * k * k + 4.0 * s1 * k, 2.0) + 16.0 * kappa * kappa * k * k)) * 0.5);
    N = floor(1.0 / h); // Number of gridpoints
    
    h = 1.0 / N; // Recalculate gridspacing
    
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
    lambdaSq = pow(c * k / h, 2);
    muSq = pow(k * kappa / (h * h), 2);
    
    kOh = (kappa * kappa) / (h * h);
    gOh = (c * c) / (h * h);
    
    // Simple (exponential) Bow Model
    a = 100; // Free parameter
    BM = sqrt(2 * a) * exp1(0.5);
    
    _Vb = -0.1; // Bowing speed
    _Fb = 80;  // Bowing force / total mass of bow;
    
    // Elasto-Plastic bow model
    
    //// the Contact Force (be with you) //////
    mus = 0.8; // static friction coeff
    mud = 0.35; // dynamic friction coeff (must be < mus!!) %EDIT: and bigger than 0
    strv = 0.08;      // "stribeck" velocity
    
    Fn = 1;    // Normal force
    
    fC = mud * Fn; // coulomb force
    fS = mus * Fn; // stiction force
    
    sig0 = 10000;                   // bristle stiffness
    sig1 = 0.1*sqrt(sig0);          // bristle damping
    sig2 = 0.4;                     // viscous friction term
    oOstrvSq = 1 / (strv * strv);   // One over strv^2
    z_ba = 0.7 * fC / sig0;         // break-away displacement (has to be < f_c/sigma_0!!)
    
    // Initialise variables for Newton Raphson
    tol = 1e-7;
    qPrev = 0;
    zPrev = 0;
    zDotPrev = 0;
    
    fp = 0;
    
    // FDS precalculations
    B1 = s0 * k;
    B2 = (2 * s1 * k) / (h * h);
    
    b1 = 2.0 / (k * k);
    b2 = (2 * s1) / (k * h * h);
    
    D = 1.0 / (1.0 + s0 * k);
    
    A1 = 2 - 2 * lambdaSq - 6 * muSq - 2 * B2;
    A2 = lambdaSq + 4 * muSq + B2;
    A3 = muSq;
    A4 = B1 - 1 + 2 * B2;
    A5 = B2;
    
    A1 *= D;
    A2 *= D;
    A3 *= D;
    A4 *= D;
    A5 *= D;
    
    E = k * k * (1 / h) * BM;
    E2 = k * k * (1 / h);
    
    reset();
    q = _Vb;
    std::cout << N << std::endl;
    K1 = -sig1 / (sig2 + 2 / k + 2 * s0);
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
    z = 0;
    zDot = 0;
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
    //    g.setColour(Colours::green);
    //    for (int i = 1; i < N; ++i)
    //    {
    //        g.fillEllipse(i * getWidth() / double(N) - 3, getHeight() / 2 - 3, 6, 6);
    //    }
    g.setColour(Colours::orange);
    for (int c = 0; c < cpIdx.size(); ++c)
    {
        g.drawEllipse(floor(cpIdx[c] * getWidth() / N - 5), floor(cy[c] - 5), 10, 10, 2);
    }
    
    g.setColour(Colours::yellow);
    g.fillEllipse(fp * getWidth() - 5, getHeight() / 2.0 - 5, 10, 10);
    
    // draw bow
    g.setColour(Colours::yellow);
    double opa = 90.0 / 100.0;
    if (opa >= 1.0)
    {
        g.setOpacity(1.0);
    }
    else
    {
        g.setOpacity(opa);
    }
    g.fillRect(floor(_bpX.load() * getWidth()), floor(_bpY.load() * getHeight()) - getHeight() / 2.0, 10, getHeight());
    g.setColour(Colour::greyLevel(0.5f).withAlpha(0.5f));
    for (double i = -12.0; i < 12.0; ++i)
    {
        double val = (1 - (pow(2.0, (i / 12.0)) - 1)) * getWidth() / 2.0;
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
    
    double Fb = _Fb.load();
    bowPos.store(clamp(_bpX.load() * N, 2, N - 3));
    int bp = floor(bowPos.load());
    bool isBowing = _isBowing;
    
    if (isBowing)
    {
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
    else if (bowModel == elastoPlastic)
    {
        excitation = E2 * (sig0 * z + sig1 * zDot + sig2 * q);
    }
    
    for (int l = 2; l < N - 2; ++l)
    {
        uNext[l] = A1 * u[l] + A2 * (u[l + 1] + u[l - 1]) - A3 * (u[l + 2] + u[l - 2]) + A4 * uPrev[l] - A5 * (uPrev[l + 1] + uPrev[l - 1]);
    }
    
    if (isBowing)
    {
        double alpha = bowPos - floor(bowPos);
        //        if (t % 10000 == 0 && stringID == 0)
        //            std::cout << alpha << " " << excitation * 100000 << std::endl;
        //        ++t;
        if (interpolation == noStringInterpol)
        {
            uNext[bp] = uNext[bp] - excitation;
        }
        else if (interpolation == linear)
        {
            uNext[bp] = uNext[bp] - excitation * (1-alpha);
            
            if (bp < N - 3)
                uNext[bp + 1] = uNext[bp + 1] - excitation * alpha;
        }
        else if (interpolation == cubic)
        {
            if (bp > 3)
                uNext[bp - 1] = uNext[bp - 1] - excitation * (alpha * (alpha - 1) * (alpha - 2)) / -6.0;
            
            uNext[bp] = uNext[bp] - excitation * ((alpha - 1) * (alpha + 1) * (alpha - 2)) / 2.0;
            
            if (bp < N - 3)
                uNext[bp + 1] = uNext[bp + 1] - excitation * (alpha * (alpha + 1) * (alpha - 2)) / -2.0;
            if (bp < N - 4)
                uNext[bp + 2] = uNext[bp + 2] - excitation * (alpha * (alpha + 1) * (alpha - 1)) / 6.0;
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
    
    b = 2.0 / k * Vb - b1 * (uI - uIPrev) - gOh * (uI1 - 2 * uI + uIM1) + kOh * (uI2 - 4 * uI1 + 6 * uI - 4 * uIM1 + uIM2) + 2 * s0 * Vb - b2 * ((uI1 - 2 * uI + uIM1) - (uIPrev1 - 2 * uIPrev + uIPrevM1));
    eps = 1;
    int i = 0;
//    std::cout << fS << std::endl;
    if (bowModel == exponential)
    {
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
//        std::cout << "Elasto-Plastic model" << std::endl;
        while (eps > tol && i < 50)
        {
            espon = exp1(-((q * q) * oOstrvSq));         //exponential function
            zss = sgn(q) * (fC + (fS - fC) * espon) / sig0;   //steady state curve: z_ss(v)
//            std::cout << zss << std::endl;
            if (q==0)
                zss = fs / sig0;
            
            // elasto-plastic function \alpha (v,z)
            alpha=0;
            if (sgn(z)==sgn(q))
            {
                if ((abs(z)>z_ba) && (abs(z)<zss))
                {
                    arg = double_Pi * (z - 0.5 * (zss + z_ba)) / (zss - z_ba);
                    alpha = 0.5 * (1 + sin(arg));
                }
                else if (abs(z)>zss)
                {
                    alpha=1;
                }
            }
            
            // non-linear function estimate
            fnl = q * (1 - alpha * z / zss);
            
            // compute derivatives
            
            // dz_ss/dv
            dz_ss = (-2 * q * sgn(q) / (strv * strv * sig0)) * (fS-fC) * espon;
            
            dalpha_v=0; //d(alpha)/dv
            dalpha_z=0; //d(alpha)/dz
            if ((sgn(z)==sgn(q)) && (abs(z)>z_ba) && (abs(z)<zss) )
            {
                double cosarg = cos(arg);
                dalpha_v = 0.5 * double_Pi * cosarg * dz_ss * (z_ba - z) / ((zss - z_ba) * (zss - z_ba));
                dalpha_z=0.5 * double_Pi * cosarg / (zss - z_ba);
            }
            
            d_fnlv = 1 - z * ((alpha + q * dalpha_v) * zss - dz_ss * alpha * q)/(zss * zss);
            d_fnlz = -q / zss * (z * dalpha_z + alpha);
            d_fnl = d_fnlv * K1 + d_fnlz * 0.5 * k;
            
            zDotNext = zDot - (fnl - zDot)/(d_fnl - 1);
            eps = abs (zDotNext-zDot);
            zDot = zDotNext;
            
            z = zPrev + k / 2 * zDotPrev + k / 2 * zDot;
            q = (-sig0 * z - sig1 * zDot - b) / (sig2 + 2/k + 2*s0);
            i = i + 1;
            if (i == 49)
                std::cout << "i = " << i << std::endl;
        }
//        std::cout << i << std::endl;
        zPrev = z;
        zDotPrev = zDot;
    }
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
    uNext = &uVecs[uNextPtrIdx][0];
    uNextPtrIdx = (uNextPtrIdx + 1) % 3;
}

Path ViolinString::generateStringPathAdvanced()
{
    auto stringBounds = getHeight() / 2.0;
    Path stringPath;
    stringPath.startNewSubPath(0, stringBounds);
    
    auto spacing = getWidth() / double(N);
    auto x = spacing;
    
    for (int y = 0; y < N; y++)
    {
        int visualScaling = (bowModel == elastoPlastic ? 500000 : 10000) * visualScale;
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
    std::cout << (bowModel == exponential ? "Exponential" : "Elasto-Plastic") << std::endl;
    if (ModifierKeys::getCurrentModifiers() == ModifierKeys::leftButtonModifier)
    {
        _Vb = 0.1;
        _Fb = 80;
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
    double maxVb = -0.2;
    
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
        float bowVelocity = e.y / (static_cast<double>(getHeight())) * maxVb;
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
    double val1 = uVec[bp - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0
    + uVec[bp] * ((alph - 1) * (alph + 1) * (alph - 2)) / 2.0
    + uVec[bp + 1] * (alph * (alph + 1) * (alph - 2)) / -2.0
    + uVec[bp + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
    
    double val = 0;
    if (bp > 3)
        val = val + uVec[bp - 1] * (alph * (alph - 1) * (alph - 2)) / -6.0;
    
    val = val + uVec[bp] * ((alph - 1) * (alph + 1) * (alph - 2)) / 2.0;
    
    if (bp < N - 3)
        val = val + uVec[bp + 1] * (alph * (alph + 1) * (alph - 2)) / -2.0;
    if (bp < N - 4)
        val = val + uVec[bp + 2] * (alph * (alph + 1) * (alph - 1)) / 6.0;
//    std::cout << val - val1 << std::endl;
    return val;
}
