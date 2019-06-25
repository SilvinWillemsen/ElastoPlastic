/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent() : minOut(-1.0), maxOut(1.0),
                                 forceSlider (Slider::RotaryVerticalDrag, Slider::TextBoxBelow),
                                 scaleGraphics (Slider::RotaryVerticalDrag, Slider::TextBoxBelow),
                                 noiseFactor (Slider::RotaryVerticalDrag, Slider::TextBoxBelow)
{
    // Make sure you set the size of the component after
    // you add any child components.
    // specify the number of input and output channels that we want to open
    setAudioChannels (0, 2);
    forceSlider.setRange (0.0, 10.0);
    forceSlider.setValue (1.0);
    addAndMakeVisible (forceSlider);
    forceSlider.addListener (this);
    
    scaleGraphics.setRange (0.01, 5.0);
    scaleGraphics.setValue (1.0);
    addAndMakeVisible (scaleGraphics);
    scaleGraphics.addListener (this);
    
    noiseFactor.setRange (0.0, 5.0);
    noiseFactor.setValue (1.0);
    addAndMakeVisible (noiseFactor);
    noiseFactor.addListener (this);
    
    noiseLabel.setText ("Noise", dontSendNotification);
    addAndMakeVisible (noiseLabel);

    forceLabel.setText ("Fn Mouse", dontSendNotification);
    addAndMakeVisible (forceLabel);
    
    scaleLabel.setText ("Scale", dontSendNotification);
    addAndMakeVisible (scaleLabel);
    
    toggleGraphics.setButtonText ("Toggle Graphics");
    addAndMakeVisible (toggleGraphics);
    toggleGraphics.addListener (this);
    toggleGraphics.setToggleState (initGraphics, sendNotification);
    
    overrideNoiseButton.setButtonText ("Override Noisecontrol");
    addAndMakeVisible (overrideNoiseButton);
    overrideNoiseButton.addListener (this);
    overrideNoiseButton.setToggleState (overrideNoise, sendNotification);
}

MainComponent::~MainComponent()
{
    HighResolutionTimer::stopTimer();
    Timer::stopTimer();
    for (auto sensel : sensels)
    {
        if (sensel->senselDetected)
        {
            sensel->shutDown();
        }
    }
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    fs = sampleRate;
    bufferSize = samplesPerBlockExpected;
    
    for (int i = 0; i < numStrings; ++i)
    {
        violinStrings.add (new ViolinString (196.0 * pow (2, (7.0 * i) / 12.0), fs, 0, elastoPlastic));
//        violinStrings.add (new ViolinString (196.0, fs, 0, elastoPlastic));
    }
//    violinStrings.add (new ViolinString (196.0, fs, 0, elastoPlastic));
//    violinStrings.add (new ViolinString (196.0, fs, 0, exponential));
//    violinStrings.add (new ViolinString (196.0, fs, 0, elastoPlastic));
    
    if (showData)
        drawData(true);
    
    setSize (appWidth, appHeight);
    
    if (toggleGraphics.getToggleStateValue() == true)
    {
        Timer::startTimerHz (15);
    }
    
    for (int i = 0; i < amountOfSensels; ++i)
    {
        sensels.add (new Sensel(0));
    }
    
    // start the hi-res timer
    if (sensels.size() != 0)
        if (sensels[0]->senselDetected)
            HighResolutionTimer::startTimer(1000.0 / 150.0);
    
    for (auto violinString : violinStrings)
    {
        addAndMakeVisible (violinString);
        std::cout << (violinString->getModel() == exponential ? "Exponential" : "ElastoPlastic") << std::endl;
    }
    if (showData)
        dataVisuals[0]->setDataPointX (violinStrings[0]->getVb());
}

void MainComponent::getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    float *const channelData1 = bufferToFill.buffer->getWritePointer(0, bufferToFill.startSample);
    float *const channelData2 = bufferToFill.buffer->getWritePointer(1, bufferToFill.startSample);
    
//    float output{0.0, 0.0};
    for (int i = 0; i < bufferSize; ++i)
    {
        float output = 0.0;
        for (auto violinString : violinStrings)
        {
            violinString->bow();
            output = output + violinString->getOutput(0.8) * (violinString->getModel() == exponential ? 800 : 80000);
            violinString->updateUVectors();
        }
        channelData1[i] = clip(output);
        channelData2[i] = clip(output);
    }
//    std::cout << channelData1[10] << std::endl;
}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
//    std::cout << violinStrings[0]->getStateAt(1) << std::endl;
    if (showData)
        drawData();
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
    
    Rectangle<int> totalArea = getLocalBounds();

    if (showControls)
    {
        Rectangle<int> controlsRect = totalArea.removeFromRight(controlsWidth);
        forceLabel.setBounds (controlsRect.removeFromTop(30));
        forceSlider.setBounds (controlsRect.removeFromTop(controlsWidth * 0.8));
        
        controlsRect.removeFromTop(20);
        scaleLabel.setBounds (controlsRect.removeFromTop(30));
        scaleGraphics.setBounds (controlsRect.removeFromTop (controlsWidth * 0.8));
        
        controlsRect.removeFromTop(20);
        noiseLabel.setBounds (controlsRect.removeFromTop(30));
        noiseFactor.setBounds (controlsRect.removeFromTop (controlsWidth * 0.8));
    
        overrideNoiseButton.setBounds (controlsRect.removeFromTop(50));
        toggleGraphics.setBounds (controlsRect.removeFromTop(50));
    }
    
    int i = 0;
    int div = appHeight / violinStrings.size();
    for (auto violinString : violinStrings)
    {
        violinString->setBounds (totalArea.removeFromTop(div));
        ++i;
    }
    
    int horSizeData = 300;
    
    if (showData)
    {
        dataVisuals[0]->setBounds (0, 0, horSizeData, 200);
        dataVisuals[1]->setBounds (horSizeData, 0, horSizeData, 200);
    }
}

void MainComponent::hiResTimerCallback()
{
    double maxVb = 0.3;
    double maxFn = 10;
    for (auto sensel : sensels)
    {
        double finger0X = 0;
        double finger0Y = 0;
        if (sensel->senselDetected)
        {
            sensel->check();
            unsigned int fingerCount = sensel->contactAmount;
            int index = sensel->senselIndex;
            for (auto violinString : violinStrings)
            {
                violinString->setBow (false);
            }
            for (int f = 0; f < fingerCount; f++)
            {
                bool state = sensel->fingers[f].state;
                float x = sensel->fingers[f].x;
                float y = sensel->fingers[f].y;
                float Vb = sensel->fingers[f].delta_y * 0.2;
                float Fn = sensel->fingers[f].force * 50;
//                std::cout << "Vb = " << Vb << " Fn = " << Fn << std::endl;
                int fingerID = sensel->fingers[f].fingerID;
                
                if (f == 0 && state) //fingerID == 0)
                {
                    finger0X = x;
                    finger0Y = y;
                    int numPlayedStrings = fingerCount;
                    if (fingerCount > violinStrings.size())
                        numPlayedStrings = violinStrings.size();
                    
                    for (int i = 0; i < numPlayedStrings; ++i)
                    {
                        Vb = clip (Vb, -maxVb, maxVb);
                        Fn = clip (Fn, 0, maxFn);
                        if (!overrideNoise)
                        {
                            violinStrings[i]->setNoise (Fn * 0.05);
                        }
//                        if (abs(Vb) == maxVb)
//                        std::cout << Vb <<std::endl;
                        violinStrings[i]->setVb (Vb);
                        violinStrings[i]->setFn (Fn);
                        violinStrings[i]->setBowPos (x, y);
                        violinStrings[i]->setBow (state);
                    }
                }
                else if (fingerID > 0 && numStrings == 1)
                {
                    //                    float dist = sqrt ((finger0X - x) * (finger0X - x) + (finger0Y - y) * (finger0Y - y));
                    float verDist = std::abs(finger0Y - y);
                    float horDist = std::abs(finger0X - x);
                    //                    std::cout << horDist << std::endl;
                    if (!(verDist <= 0.3 && horDist < 0.05))
                    {
                        violinStrings[index]->setFingerPosition(x);
                    }
                }
            }
        }
    }
}

void MainComponent::timerCallback()
{
    repaint();
}

double MainComponent::clip (double output, double min, double max)
{
    if (output > max)
    {
        return output = max;
    }
    else if (output < min)
    {
        return output = min;
    }
    return output;
}

void MainComponent::sliderValueChanged (Slider* slider)
{
    if (slider == &forceSlider)
    {
        for (auto violinString : violinStrings)
        {
            violinString->setFn (slider->getValue());
        }
    }
    
    if (slider == &scaleGraphics)
    {
        for (auto violinString : violinStrings)
        {
            violinString->scaleVisuals (slider->getValue());
        }
    }
    
    if (slider == &noiseFactor)
    {
        for (auto violinString : violinStrings)
        {
            violinString->setNoise (slider->getValue());
        }
    }
}

void MainComponent::buttonClicked (Button* button)
{
    if (button == &toggleGraphics)
    {
        if (toggleGraphics.getToggleStateValue() == false)
            Timer::stopTimer();
        else
            Timer::startTimerHz(60);
    }
    
    if (button == &overrideNoiseButton)
    {
        overrideNoise = button->getToggleState();
        if (overrideNoise)
        {
            sliderValueChanged (&noiseFactor);
        }
    }
}

void MainComponent::drawData(bool createNew)
{
    float detail = 1000.0;
    
    std::vector<double> data (detail, 0);
    std::vector<double> xData (detail, 0);
    for (int i = 0; i < detail; ++i)
    {
        double q = (i - detail * 0.5) / detail;
        double espon = exp(-((q * q) / (0.1 * 0.1)));
        data[i] = sgn(q) * (violinStrings[0]->getFC() + (violinStrings[0]->getFS() - violinStrings[0]->getFC()) * espon) / 10000.0;
        xData[i] = q;
    }
    if (createNew)
    {
        dataVisuals.add (new DataVisual (steadyState, data, xData, 0, true, true));
        addAndMakeVisible (dataVisuals[0]);
    }
    dataVisuals[0]->setData (data);
//    dataVisuals[0]->setXData (xDataAlpha);
    dataVisuals[0]->setDataPointX (violinStrings[0]->getVb());
    
    std::vector<double> dataAlpha (detail, 0);
    std::vector<double> xDataAlpha (detail, 0);
    double zssPlot = abs(data[dataVisuals[0]->getVbIdx()]);
    double z_ba = 0.7 * violinStrings[0]->getFC() / violinStrings[0]->getSig0();
    int j = 0;
    for (int i = 0; i < detail; ++i)
    {
        double zVec = (i - detail * 0.5) / (detail * 100.0);
        xDataAlpha[i] = zVec;
        if ((abs(zVec)>z_ba) && (abs(zVec)<zssPlot))
        {
            ++j;
            dataAlpha[i]=0.5*(1+sin(sgn(zVec)*(double_Pi*(zVec-sgn(zVec)*0.5*(zssPlot+z_ba))/(zssPlot-z_ba))));
        }
        else if (abs(zVec)>zssPlot)
            dataAlpha[i]=1;
    }
    if (createNew)
    {
        dataVisuals.add (new DataVisual (adhesionMap, dataAlpha, xDataAlpha, 0, true, true));
        addAndMakeVisible (dataVisuals[1]);
    }
//    if (j == 0)
//    {
//        std::cout << "wait" << std::endl;
//    }
    dataVisuals[1]->setData (dataAlpha);
    dataVisuals[1]->setXData (xDataAlpha);
    dataVisuals[1]->setDataPointX(violinStrings[0]->getZ());
    double curZ = violinStrings[0]->getZ();
    double valueForZ = 0;
    if ((abs(curZ)>z_ba) && (abs(curZ)<zssPlot))
    {
        valueForZ=0.5*(1+sin(sgn(curZ)*(double_Pi*(curZ-sgn(curZ)*0.5*(zssPlot+z_ba))/(zssPlot-z_ba))));
    }
    else if (abs(curZ)>zssPlot)
        valueForZ = 1;
    dataVisuals[1]->setDataPointY (valueForZ);
}
