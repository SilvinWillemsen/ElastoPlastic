/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "ViolinString.h"
#include "../SenselWrapper/SenselWrapper.h"
//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent   : public AudioAppComponent,
                        public HighResolutionTimer,
                        public Timer
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent();

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (Graphics& g) override;
    void resized() override;
    
    void hiResTimerCallback() override;
    void timerCallback() override;
    
    float clip(float output);
    
//    void sliderValueChanged (Slider* slider) override;
//    void buttonClicked (Button* button) override;

private:
    //==============================================================================
    // Your private member variables go here...
    double fs;
    double bufferSize;
    
    float minOut;
    float maxOut;
    
    int numStrings;
    OwnedArray<ViolinString> violinStrings;
    OwnedArray<Sensel> sensels;
    
    unsigned long stateUpdateCounter = 0;
    
    int appWidth = 1440;
    int appHeight = 800;
    int controlsWidth = 100;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
