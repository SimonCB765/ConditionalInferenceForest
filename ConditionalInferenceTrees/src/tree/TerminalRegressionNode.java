/**
 * 
 */
package tree;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Simon
 *
 */
public class TerminalRegressionNode extends Node
{

	double meanResponseValue;

	TerminalRegressionNode(String loadString)
	{
		loadString = loadString.replaceAll("\n", "");
		String split[] = loadString.split("\t");
		this.nodeDepth = Integer.parseInt(split[0]);
		this.numberOfObservationsInNode = Integer.parseInt(split[1]);
		this.meanResponseValue = Double.parseDouble(split[2]);
		this.pValue = 0.0;
		pValues = new HashMap<String, Double>();
		testStatistics = new HashMap<String, Double>();
		pValues.put("NA", 0.0);
		testStatistics.put("NA", 0.0);
	}

	TerminalRegressionNode(int currentDepth, int numberOfObservationsInNode, double meanResponseValue)
	{
		this.nodeDepth = currentDepth;
		this.numberOfObservationsInNode = numberOfObservationsInNode;
		this.meanResponseValue = meanResponseValue;
		this.pValue = 0.0;
		pValues = new HashMap<String, Double>();
		testStatistics = new HashMap<String, Double>();
		pValues.put("NA", 0.0);
		testStatistics.put("NA", 0.0);
	}

	void display()
	{
		for (int i = 0; i < nodeDepth; i++)
		{
			System.out.print("|  ");
		}
		System.out.format("Mean value : %f, Observations : %d.\n", meanResponseValue, numberOfObservationsInNode);
	}

	Object predict(Map<String, Object> currentObservation)
	{
		return predict(currentObservation, false);
	}

	Object predict(Map<String, Object> currentObservation, boolean isProbabilisticPrediction)
	{
		return meanResponseValue;
	}

	String save()
	{
		String returnValue = Integer.toString(this.nodeDepth) + "\t" + Integer.toString(this.numberOfObservationsInNode) + "\t" +
				Double.toString(meanResponseValue);
		return returnValue;
	}

}
