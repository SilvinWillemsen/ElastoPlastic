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
    
    g.fillAll (Colours::black);
    if (isResized)
    {
        if (drawGraph)
        {
            for (int i = 0; i < numDataPoints; ++i)
            {
                normData[i] = height - height * (data[i] - minVal) / (maxVal - minVal) + halfMargin;
                normXData[i] = width * (xData[i] - minXVal) / (maxXVal - minXVal) + halfMargin;
            }
//            g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));   // clear the background
            g.setColour(Colours::blue);
            Path graph;
            
            graph.startNewSubPath(normXData[0], normData[0]);
            for (int i = 1; i < numDataPoints; i++)
            {
                graph.lineTo(normXData[i], normData[i]);
//                if (data[i] != 1 && data[i] != 0)
//                {
//                    std::cout << "i = " << i << " data = " << data[i] << std::endl;
//                }
            }
//            graph.lineTo(width - halfMargin, normData[numDataPoints - 1]);
            g.strokePath(graph, PathStrokeType(2.0f));
        }
        if (drawDataPoint)
        {
            g.setColour(Colours::orange);
            std::ostringstream dpValX;
            std::ostringstream dpValY;
            if (graphType == steadyState)
            {
                int xPos = width * (dataPointX - minXVal) / (maxXVal - minXVal) + halfMargin;
                auto testIt = std::find(xData.begin(), xData.end(), round(dataPointX * numDataPoints) / static_cast<float> (numDataPoints));
                int idx = std::distance(xData.begin(), testIt);
                g.drawEllipse (xPos - 5, normData[idx] - 5, 10, 10, 2);
                dpValX << std::setprecision (2) << dataPointX;
                dpValY << std::setprecision (2) << data[idx];
                dataPointLabel.setText("v = " + dpValX.str()+ "\n" + "f = " + dpValY.str(), dontSendNotification);
            }
            else if (graphType == adhesionMap)
            {
                int xPos = width * (dataPointX - minXVal) / (maxXVal - minXVal) + halfMargin;
                int yPos = height * (dataPointY - minVal) / (maxVal - minVal) - halfMargin;
                g.drawEllipse (xPos - 5, height - yPos - 5, 10, 10, 2);
                dpValX << std::setprecision (2) << dataPointX;
                dpValY << std::setprecision (2) << dataPointY;
                juce::String alphaString = String(CharPointer_UTF8 ("\xce\xb1")) + " = " + dpValY.str() + "\n" + "z = " + dpValX.str();
                juce::String totString = juce::String(alphaString) ;
                dataPointLabel.setText(totString, dontSendNotification);
            }

        }
    }
    g.setColour (Colours::white);
    g.drawRect (getLocalBounds(), 2);
}

void DataVisual::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..
    width = getWidth() - margin;
    height = getHeight() - margin;

    for (int i = 0; i < numDataPoints; ++i)
    {
        normData[i] = height - height * (data[i] - minVal) / (maxVal - minVal) + halfMargin;
        normXData[i] = width * (xData[i] - minXVal) / (maxXVal - minXVal) + halfMargin;
    }
    
    if (drawDataPoint)
    {
        dataPointLabel.setBounds(getWidth() - margin - 100, getHeight() - margin - 30, 100, 30);
        dataPointLabel.setJustificationType (Justification::right);
        addAndMakeVisible (dataPointLabel);
    }
    isResized = true;
}
