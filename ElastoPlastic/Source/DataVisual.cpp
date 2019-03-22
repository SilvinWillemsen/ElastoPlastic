/*
  ==============================================================================

    DataVisual.cpp
    Created: 21 Mar 2019 2:11:24pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include "../JuceLibraryCode/JuceHeader.h"
#include "DataVisual.h"

//==============================================================================
DataVisual::DataVisual(GraphType graphType, std::vector<double> data, std::vector<double> xData,
                       double dataPointX, bool drawGraph, bool drawDataPoint) : graphType (graphType),
                                                                                data (data),
                                                                                xData (xData),
                                                                                dataPointX (dataPointX),
                                                                                drawGraph (drawGraph),
                                                                                drawDataPoint (drawDataPoint),
                                                                                numDataPoints (data.size())
{
    minVal = *std::min_element (data.begin(), data.end());
    maxVal = *std::max_element(data.begin(), data.end());
    minXVal = *std::min_element(xData.begin(), xData.end());
    maxXVal = *std::max_element(xData.begin(), xData.end());
    normData.resize (numDataPoints);
    normXData.resize (numDataPoints);
}

DataVisual::~DataVisual()
{
}

void DataVisual::paint (Graphics& g)
{
    /* This demo code just fills the component's background and
       draws some placeholder text to get you started.

       You should replace everything in this method with your own
       drawing code..
    */
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    if (isResized)
    {
        if (drawGraph)
        {
            for (int i = 0; i < numDataPoints; ++i)
            {
                normData[i] = height - height * (data[i] - minVal) / (maxVal - minVal);
                normXData[i] = width * (xData[i] - minXVal) / (maxXVal - minXVal);
            }
            g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));   // clear the background
            g.setColour(Colours::yellow);
            Path graph;
            
            graph.startNewSubPath(0, normData[0]);
            for (int i = 1; i < numDataPoints - 1; i++)
            {
                graph.lineTo(normXData[i], normData[i]);
//                if (data[i] != 1 && data[i] != 0)
//                {
//                    std::cout << "i = " << i << " data = " << data[i] << std::endl;
//                }
            }
            graph.lineTo(width, normData[numDataPoints - 1]);
            g.strokePath(graph, PathStrokeType(2.0f));
        }
        if (drawDataPoint)
        {
            g.setColour(Colours::orange);
            if (graphType == steadyState)
            {
                int xPos = width * (dataPointX - minXVal) / (maxXVal - minXVal);
                auto testIt = std::find(xData.begin(), xData.end(), round(dataPointX * numDataPoints) / static_cast<float> (numDataPoints));
                int idx = std::distance(xData.begin(), testIt);
                g.drawEllipse (xPos - 5, normData[idx] - 5, 10, 10, 2);
            } else if (graphType == adhesionMap)
            {
                int xPos = width * (dataPointX - minXVal) / (maxXVal - minXVal);
                int yPos = height * (dataPointY - minVal) / (maxVal - minVal);
                g.drawEllipse (xPos - 5, height - (yPos) - 5, 10, 10, 2); }
            
        }
    }
}

void DataVisual::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..
    width = getWidth();
    height = getHeight();

    for (int i = 0; i < numDataPoints; ++i)
    {
        normData[i] = height - height * (data[i] - minVal) / (maxVal - minVal);
        normXData[i] = width * (xData[i] - minXVal) / (maxXVal - minXVal);
    }
    
    isResized = true;
}
