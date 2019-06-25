/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "ViolinString.h"
#include "../SenselWrapper/SenselWrapper.h"
#include "DataVisual.h"
//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent   : public AudioAppComponent,
                        public HighResolutionTimer,
                        public Timer,
                        public Slider::Listener,
                        public Button::Listener
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
    
    double clip (double output, double min = -1.0, double max = 1.0);
    
    void sliderValueChanged (Slider* slider) override;
    void buttonClicked (Button* button) override;

    int sgn (double val) { return (0 < val) - (val < 0); };
    
    void drawData (bool createNew = false);
private:
    //==============================================================================
    // Your private member variables go here...
    double fs;
    double bufferSize;
    
    float minOut;
    float maxOut;
    
    OwnedArray<ViolinString> violinStrings;
    OwnedArray<Sensel> sensels;
    OwnedArray<DataVisual> dataVisuals;
    unsigned long stateUpdateCounter = 0;
    
    int appWidth = 1440;
    int appHeight = 800;
    int controlsWidth = 100;

    Slider forceSlider;
    Slider scaleGraphics;
    Slider noiseFactor;
    
    ToggleButton toggleGraphics;
    ToggleButton overrideNoiseButton;
    Label forceLabel;
    Label scaleLabel;
    Label noiseLabel;
    bool initGraphics = true;
    
    int amountOfSensels = 1;
    
    bool showData = true;
    bool showControls = true;
    int numStrings = 4;
    bool overrideNoise = false;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
