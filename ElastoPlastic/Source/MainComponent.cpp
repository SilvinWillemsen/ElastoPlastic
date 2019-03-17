/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent() : minOut(-1.0), maxOut(1.0)
{
    // Make sure you set the size of the component after
    // you add any child components.
    // specify the number of input and output channels that we want to open
    setAudioChannels (0, 1);
}

MainComponent::~MainComponent()
{
//    for (auto sensel : sensels)
//    {
//        if (sensel->senselDetected)
//        {
//            sensel->shutDown();
//        }
//    }
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    fs = sampleRate;
    bufferSize = samplesPerBlockExpected;
    
    violinStrings.add (new ViolinString (176.0, fs, 0, elastoPlastic));
    violinStrings.add (new ViolinString (196.0, fs, 0, exponential));
    setSize (appWidth, appHeight);
    Timer::startTimerHz(60);
    
    for (auto violinString : violinStrings)
    {
        addAndMakeVisible (violinString);
    }
}

void MainComponent::getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    float *const channelData1 = bufferToFill.buffer->getWritePointer(0, bufferToFill.startSample);
//    float *const channelData2 = bufferToFill.buffer->getWritePointer(1, bufferToFill.startSample);
    
//    float output{0.0, 0.0};
    for (int i = 0; i < bufferSize; ++i)
    {
        float output = 0.0;
        for (auto violinString : violinStrings)
        {
            violinString->bow();
            output = output + violinString->getOutput(0.3) * (violinString->getModel() == exponential ? 800 : 3000);
            violinString->updateUVectors();
        }
        channelData1[i] = clip(output);
//        channelData2[i] = clip(output);
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
    int i = 0;
    int div = appHeight / violinStrings.size();
    for (auto violinString : violinStrings)
    {
        violinString->setBounds (0, div * i, appWidth, div);
        ++i;
    }
}

void MainComponent::hiResTimerCallback()
{
    //do senselstuff
}

void MainComponent::timerCallback()
{
    repaint();
}

float MainComponent::clip(float output)
{
    if (output > maxOut)
    {
        return output = maxOut;
    }
    else if (output < minOut)
    {
        return output = minOut;
    }
    return output;
}
