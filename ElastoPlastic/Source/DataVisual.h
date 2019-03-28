/*
  ==============================================================================

    DataVisual.h
    Created: 21 Mar 2019 2:11:24pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"

//==============================================================================
/*
*/
enum GraphType
{
    steadyState,
    adhesionMap,
};

class DataVisual    : public Component
{
public:
    DataVisual(GraphType graphType, std::vector<double> data, std::vector<double> xData, double dataPointX, bool drawGraph = true, bool drawDataPoint = false);
    ~DataVisual();

    void paint (Graphics&) override;
    void resized() override;

    void setDataPointX (double dpX) { dataPointX = dpX; };
    void setDataPointY (double dpY) { dataPointY = dpY; };
    void setData (std::vector<double> d) { data = d; };
    void setXData (std::vector<double> xd) { xData = xd; };
    
    int getVbIdx() { return VbIdx; };
private:
    std::vector<double> data;
    std::vector<int> normData;
    std::vector<double> xData;
    std::vector<int> normXData;

    double dataPointX;
    double dataPointY;
    int normDataPointX;
    int normDataPointY;
    int width = 100;
    int height = 100;
    double minVal;
    double maxVal;
    double minXVal;
    double maxXVal;
    bool isResized = false;
    
    bool drawGraph;
    bool drawDataPoint;
    
    int VbIdx = 0;
    unsigned long numDataPoints;
    
    GraphType graphType;
    static const int margin = 20;
    static const int halfMargin = margin * 0.5;
    
    Label dataPointLabel;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (DataVisual)
};
