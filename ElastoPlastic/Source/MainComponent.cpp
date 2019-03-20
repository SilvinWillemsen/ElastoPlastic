/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent() : minOut(-1.0), maxOut(1.0),
                                 forceSlider (Slider::RotaryVerticalDrag, Slider::TextBoxBelow),
                                 scaleGraphics (Slider::RotaryVerticalDrag, Slider::TextBoxBelow)
{
    // Make sure you set the size of the component after
    // you add any child components.
    // specify the number of input and output channels that we want to open
    setAudioChannels (0, 2);
    forceSlider.setRange (0.0, 100.0);
    forceSlider.setValue (1.0);
    addAndMakeVisible(forceSlider);
    forceSlider.addListener(this);
    forceSlider.setSliderStyle (Slider::RotaryVerticalDrag);
    
    scaleGraphics.setRange (1.0, 10.0);
    scaleGraphics.setValue (5.0);
    addAndMakeVisible(scaleGraphics);
    scaleGraphics.addListener(this);
    scaleGraphics.setSliderStyle (Slider::RotaryVerticalDrag);
    
    forceLabel.setText ("Fn Mouse", dontSendNotification);
    addAndMakeVisible (forceLabel);
    scaleLabel.setText ("Scale", dontSendNotification);
    addAndMakeVisible (scaleLabel);
    toggleGraphics.setButtonText("Toggle Graphics");
    addAndMakeVisible(toggleGraphics);
    toggleGraphics.addListener(this);
    toggleGraphics.setToggleState(initGraphics, sendNotification);
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
    
    violinStrings.add (new ViolinString (110.0, fs, 0, elastoPlastic));
//    violinStrings.add (new ViolinString (196.0, fs, 0, exponential));
    setSize (appWidth, appHeight);
    
    if (toggleGraphics.getToggleStateValue() == true)
    {
        Timer::startTimerHz (60);
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
            output = output + violinString->getOutput(0.8) * (violinString->getModel() == exponential ? 800 : 3500);
            violinString->updateUVectors();
        }
        channelData1[i] = clip(output);
        channelData2[i] = clip(output);
    }
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

    // You can add your drawing code here!
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
    
    Rectangle<int> totalArea = getLocalBounds();
    Rectangle<int> controlsRect = totalArea.removeFromRight(controlsWidth);
    
    int i = 0;
    int div = appHeight / violinStrings.size();
    for (auto violinString : violinStrings)
    {
        violinString->setBounds (totalArea.removeFromTop(div));
        ++i;
    }
    forceLabel.setBounds(controlsRect.removeFromTop(30));
    forceSlider.setBounds (controlsRect.removeFromTop(controlsWidth * 0.8));
    controlsRect.removeFromTop(20);
    scaleLabel.setBounds(controlsRect.removeFromTop(30));
    scaleGraphics.setBounds (controlsRect.removeFromTop (controlsWidth * 0.8));
    toggleGraphics.setBounds (controlsRect.removeFromTop(100));
    
}

void MainComponent::hiResTimerCallback()
{
    double maxVb = 0.2;
    double maxFn = 1000;
    for (auto sensel : sensels)
    {
        double finger0X = 0;
        double finger0Y = 0;
        if (sensel->senselDetected)
        {
            sensel->check();
            unsigned int fingerCount = sensel->contactAmount;
            int index = sensel->senselIndex;
            for (int f = 0; f < fingerCount; f++)
            {
                bool state = sensel->fingers[f].state;
                float x = sensel->fingers[f].x;
                float y = sensel->fingers[f].y;
                float Vb = sensel->fingers[f].delta_y / 10.0;
                float Fn = sensel->fingers[f].force * 1000;
                std::cout << "Vb = " << Vb << " Fn = " << Fn << std::endl;
                int fingerID = sensel->fingers[f].fingerID;
                
                if (f == 0 && state) //fingerID == 0)
                {
                    finger0X = x;
                    finger0Y = y;
                    violinStrings[index]->setBow(state);
                    Vb = clip (Vb, -maxVb, maxVb);
                    Fn = clip (Fn, 0, maxFn);
                    violinStrings[index]->setVb(Vb);
                    violinStrings[index]->setFn(Fn);
                    violinStrings[index]->setBowPos(x, y);
                }
//                else if (fingerID > 0)
//                {
//                    //                    float dist = sqrt ((finger0X - x) * (finger0X - x) + (finger0Y - y) * (finger0Y - y));
//                    float verDist = std::abs(finger0Y - y);
//                    float horDist = std::abs(finger0X - x);
//                    //                    std::cout << horDist << std::endl;
//                    if (!(verDist <= 0.3 && horDist < 0.05))
//                    {
//                        violinStrings[index]->setFingerPosition(x);
//                    }
//                }
            }
            
            if (fingerCount == 0)
            {
                violinStrings[index]->setBow(false);
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
}
