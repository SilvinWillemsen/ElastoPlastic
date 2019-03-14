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
    setAudioChannels (2, 2);
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
    
    violinStrings.add (new ViolinString (110.0, fs, 0));
    setSize (appWidth, appHeight);
    Timer::startTimerHz(15);
    
    addAndMakeVisible (violinStrings[0]);
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
    float output = 0;
    for (int i = 0; i < bufferSize; ++i)
    {
        violinStrings[0]->bow();
        output = violinStrings[0]->getOutput (0.3) * 800;
        violinStrings[0]->updateUVectors();
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
    violinStrings[0]->setBounds (0, 0, appWidth, appHeight);
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
